#!/usr/bin/env python3

# Written by Emmy Killett, Claude Opus 4.6 Extended Thinking, and Claude Opus 4.7 Adaptive.

from __future__ import annotations

import argparse
import csv
import datetime as dt
import logging
import os
import re
import shutil
import subprocess
import sys
import timeit
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, BinaryIO, Final

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import dask
import h5py  # noqa: F401 — needed by xarray for NetCDF4 backend
import matplotlib.pyplot as plt
import nfft
import numpy as np
import xarray as xr
from colour import SpectralDistribution, XYZ_to_sRGB, sd_to_XYZ
from colour.colorimetry import MSDS_CMFS_STANDARD_OBSERVER
from tqdm import tqdm

dask.config.set(scheduler="synchronous")
# Cap xarray's open-file cache.  With thousands of AQUA_MODIS L3m files,
# the default of 128 + repeated tile reads thrashes HDF5/libnetcdf metadata.
xr.set_options(file_cache_maxsize=32)
import gc

__version__ = "0.2.1"

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------

# Per-input-choice dispatch: (y_key, folder_attr_on_options, glob_pattern).
# y_key stays hardcoded because it's the science variable; the lat/lon/time
# coord names are auto-detected from the first matching file.
_INPUT_DISPATCH: Final[dict[str, tuple[str, str, str]]] = {
    "SSHA": ("SLA", "ssha_folder", "*.nc"),
    "simple_grids": ("ssha", "simple_folder", "**/*.nc"),
    "Argo": ("ohc_2d_lt_700m", "argo_folder", "*.nc"),
    "MUR_SST": ("sst_anomaly", "mur_sst_folder", "*.nc"),
    "AQUA_MODIS": ("chlor_a", "aqua_modis_folder", "*.nc"),
    "GRACE": ("lwe_thickness", "grace_folder", "*.nc"),
}


class Options:
    """All global options in one place."""

    def __init__(self) -> None:
        """Initialize Options with default values."""
        # Identity
        self.my_name: str = Path(sys.argv[0]).stem
        self.cwd: Path = Path.cwd().expanduser().resolve(strict=True)
        self.log_mode: int = logging.INFO
        self.args: argparse.Namespace = argparse.Namespace()

        # Paths
        self.output_base: Path = self.cwd / "output"
        # self.ssha_folder: Path = (
        #     self.cwd / "sealevel_spectra" / "fast_202306" / "fast_netCDF4"
        # )
        self.ssha_folder: Path = self.cwd / "sealevel_spectra" / "SSHA_v2205"
        self.simple_folder: Path = self.cwd / "sealevel_spectra" / "simple_grids"
        self.argo_folder: Path = self.cwd / "sealevel_spectra" / "Argo"
        self.mur_sst_folder: Path = (
            self.cwd / "sealevel_spectra" / "MUR_SST" / "MUR25-JPL-L4-GLOB-v04.2"
        )
        self.aqua_modis_folder: Path = self.cwd / "sealevel_spectra" / "AQUA_MODIS"
        self.grace_folder: Path = self.cwd / "sealevel_spectra" / "JPL_GRACE_mascons"
        self.cie_file: Path = (
            self.cwd / "sealevel_spectra" / "ciexyz31_1_trimmed_380nm_760nm.csv"
        )
        self.output_folder: Path = Path()  # set in main()

        # Configuration

        self.input_choices_list: list[str] = [
            "SSHA",
            "simple_grids",
            "Argo",
            "MUR_SST",
            "AQUA_MODIS",
            "GRACE",
        ]
        self.input_choice: str = "AQUA_MODIS"
        self.min_period: float = 3 * 7.0  # days
        self.max_period: float = 24 * 7.0  # days
        # self.min_period: float = 10 * 7.0  # days
        # self.max_period: float = 54 * 7.0  # days

        self.fit_terms: list[str] = []
        self.fit_terms += ["constant", "trend", "accel"]
        # self.fit_terms += ["annual"]
        # self.fit_terms += ["semiannual"]
        # self.fit_terms += ["annual", "semiannual"]

        self.xskip: int = 100  # only use every "xskip" point in each direction.
        self.tskip: int = 1  # only use every "tskip" point in the time series

        # Desired memory footprints; chunk_size and tile_side are derived from these
        self.chunk_bytes: float = 250e6  # bytes per dask time chunk
        self.block_bytes: float = 500e6  # bytes per loop (lat, lon) tile

        # Controls xarray.open_mfdataset's parallel setting.
        # Be careful, parallel=True can spike RAM on open
        # (Above, not unrelated: dask.config.set(scheduler="synchronous"))
        self.mf_parallel: bool = bool(0)

        self.thepower: float = 0.8

        self.dpi: int = 300
        self.figsize: tuple[int, int] = (10, 5)

        # Function selection
        self.use_new_funcs: bool = bool(0)

        # Files to copy into output
        self.files_to_copy: list[str] = ["projections.sh", "overflow.sh", "notation.sh"]

        # Not a configurable option.
        self.fake_chunk: int = 100

        # Hard cap on total dask chunks across all opened files.  Sets the
        # upper bound on the HighLevelGraph size at ~300 bytes/task ⇒
        # ~300 MB graph at the default. Lower if your machine is tight on
        # RAM at open time.
        self.max_graph_chunks: int = 1_000_000

        # Spatial chunks for xr.open_mfdataset; populated by
        # derive_spatial_chunks(options) before load_input_data().
        self.lat_chunk: int = -1
        self.lon_chunk: int = -1
        self.tile_side: int = 0

        # First-file shape, populated by configure_keys_for_input().
        self.n_lat_raw: int = 0
        self.n_lon_raw: int = 0
        self.on_disk_chunks: tuple[int, int] | None = None

        # Kerchunk virtual-zarr loader (opt-in).  When use_kerchunk is True,
        # multi-file datasets are read through an fsspec reference filesystem
        # over a sidecar JSON, which lets strided isel translate to chunk-
        # aware reads (real I/O reduction when xskip >= disk-chunk-size).
        self.use_kerchunk: bool = False
        self.kerchunk_path: Path | None = None

        # Data structures (initialized by populate_results after hot path)
        self.plot_options: PlotOptions = PlotOptions()
        self.grid: GridData = GridData()
        self.results: Results = Results()

        # Data key names (set per input_choice in configure_keys_for_input)
        self.time_key: str = ""
        self.y_key: str = ""
        self.lat_key: str = ""
        self.lon_key: str = ""


# ---------------------------------------------------------------------------
# Data structures (mirrors C++ structs from abridged-definitions.hpp)
# ---------------------------------------------------------------------------


@dataclass
class LatLonData:
    """Equal-grid latlon data (mirrors C++ latlon_s, definitions.hpp:238-256).

    Holds coordinates and output data for a 2D lat/lon grid.
    """

    lat: list[float] = field(default_factory=list)
    lon: list[float] = field(default_factory=list)
    outputs: list[list[list[float]]] = field(default_factory=list)  # [param][lat][lon]
    elev: list[list[float]] = field(default_factory=list)
    mascon_lats: list[float] = field(default_factory=list)
    mascon_lons: list[float] = field(default_factory=list)
    mask: float = float("nan")


@dataclass
class AnalysisOptions:
    """Analysis options (subset of C++ analysis_options_s used by GMT path).

    Controls which output format and parameter types are used.
    """

    output_choice: int = 5  # 1=unstructured, 5=latlon grid
    output_type: list[int] = field(
        default_factory=list
    )  # per-parameter (e.g. 104=phase)


@dataclass
class PlotOptions:
    """Plotting/visualization options (mirrors C++ plot_options_s, definitions.hpp:443-490).

    Controls GMT script generation, projection, color scheme, and matplotlib output.
    """

    outputfolder: str = ""
    just_the_filenames: list[str] = field(default_factory=list)
    output_base: str = "map_parameter"
    projection: int = 2  # 2=Winkel Tripel (matches C++ default)
    coastlines: int = 1  # 1=coast, 2=coast+InSAR
    symmetric_limit: float = -1.0  # <=0 disables (matches C++)
    scale_digits: int = 2
    color_scheme: int = 2  # 1=white BG, 2=black BG
    montage: int = 0
    blurb_disabled: int = 1
    phase_mask: int = 12  # 1/2=abrupt, 11/12=gradual (matches C++)
    plot_mascons: int = 0
    no_gmt_plots: int = 1  # 1=skip GMT
    # Matplotlib-only fields:
    dpi: int = 300
    show_fig: bool = False
    save_fig: bool = True
    figsize: tuple[int, int] = (14, 7)
    # Disabled land and ocean fill because it only seems to happen near the poles.
    # land_color: str = "#333333"  # Dark grey
    # ocean_color: str = "black"  # Matches dark_background theme
    dark_mode: bool = bool(1)  # False = light/white background
    show_borders: bool = False  # Political boundaries
    show_rivers: bool = False  # Major rivers
    coastline_resolution: str = "110m"  # "110m" or "50m"


@dataclass
class GridData:
    """Grid point data (simplified from C++ grid_s, definitions.hpp:258-319).

    Holds grid coordinates and boundaries.
    """

    lat: list[float] = field(default_factory=list)
    lon: list[float] = field(default_factory=list)
    maxlat: float = 90.0
    minlat: float = -90.0
    maxlon: float = 360.0
    minlon: float = 0.0
    latlon: LatLonData | None = None


@dataclass
class Results:
    """Computation results (mirrors C++ results_s, definitions.hpp:338-441).

    Holds output data, metadata, spatial extent, and RGB color mapping state.
    """

    titles: list[str] = field(default_factory=list)
    units: list[str] = field(default_factory=list)
    base_unit: str = "cm"
    latlon: LatLonData = field(default_factory=LatLonData)
    error_bars: list[float] = field(default_factory=list)
    marker_lats: list[list[float]] = field(default_factory=list)
    marker_lons: list[list[float]] = field(default_factory=list)
    minlat: float = -90.0
    maxlat: float = 90.0
    minlon: float = 0.0
    maxlon: float = 360.0
    rgb: list[float] = field(default_factory=list)
    rgb_choice: int = 2  # C++ default: 0=non-uniform, 1=uniform, 2=multi-width
    max_labels: int = 10  # C++ default
    max_widths: int = 1
    options: AnalysisOptions = field(default_factory=AnalysisOptions)


DAYS_IN_YEAR: Final[float] = 365.25

# XYZ-to-linear-sRGB matrix (Hughes & Williams 2010, no gamma correction)
A_OLD_XYZ_TO_SRGB: Final[np.ndarray] = np.array(
    [
        [3.2409699, -1.5373832, -0.49861079],
        [-0.96924375, 1.8759676, 0.041555082],
        [0.055630032, -0.20397685, 1.0569714],
    ]
)

# ---------------------------------------------------------------------------
# Coordinate auto-detection
# ---------------------------------------------------------------------------

_LATITUDE_NAMES: Final[frozenset[str]] = frozenset(
    {"lat", "latitude", "Latitude", "LATITUDE", "nav_lat", "y"}
)
_LONGITUDE_NAMES: Final[frozenset[str]] = frozenset(
    {"lon", "longitude", "Longitude", "LONGITUDE", "nav_lon", "x"}
)
_TIME_NAMES: Final[frozenset[str]] = frozenset({"time", "Time", "TIME", "t"})

_LATITUDE_UNITS: Final[frozenset[str]] = frozenset(
    {"degrees_north", "degree_north", "degrees_N", "degree_N"}
)
_LONGITUDE_UNITS: Final[frozenset[str]] = frozenset(
    {"degrees_east", "degree_east", "degrees_E", "degree_E"}
)
_TIME_SINCE_PATTERN: Final[re.Pattern[str]] = re.compile(r"\bsince\b", re.IGNORECASE)


def _detect_latitude(ds: xr.Dataset) -> str | None:
    """Detect the latitude variable in an xarray Dataset.

    Uses a CF-style priority cascade: standard_name > units > conventional
    name > axis attribute.

    Args:
        ds: An open xarray Dataset.

    Returns:
        The name of the latitude variable, or None if no match is found.
    """
    for name, var in ds.variables.items():
        if var.attrs.get("standard_name") == "latitude":
            return str(name)
    for name, var in ds.variables.items():
        if str(var.attrs.get("units", "")) in _LATITUDE_UNITS:
            return str(name)
    for name in ds.variables:
        if str(name) in _LATITUDE_NAMES:
            return str(name)
    for name, var in ds.variables.items():
        if var.attrs.get("axis") == "Y":
            return str(name)
    return None


def _detect_longitude(ds: xr.Dataset) -> str | None:
    """Detect the longitude variable in an xarray Dataset.

    Uses a CF-style priority cascade: standard_name > units > conventional
    name > axis attribute.

    Args:
        ds: An open xarray Dataset.

    Returns:
        The name of the longitude variable, or None if no match is found.
    """
    for name, var in ds.variables.items():
        if var.attrs.get("standard_name") == "longitude":
            return str(name)
    for name, var in ds.variables.items():
        if str(var.attrs.get("units", "")) in _LONGITUDE_UNITS:
            return str(name)
    for name in ds.variables:
        if str(name) in _LONGITUDE_NAMES:
            return str(name)
    for name, var in ds.variables.items():
        if var.attrs.get("axis") == "X":
            return str(name)
    return None


def _detect_time(ds: xr.Dataset) -> str | None:
    """Detect the time variable in an xarray Dataset.

    Uses a CF-style priority cascade: standard_name > "since" units pattern >
    conventional name > axis attribute.

    Args:
        ds: An open xarray Dataset.

    Returns:
        The name of the time variable, or None if no match is found.
    """
    for name, var in ds.variables.items():
        if var.attrs.get("standard_name") == "time":
            return str(name)
    for name, var in ds.variables.items():
        units = var.attrs.get("units")
        if units is not None and _TIME_SINCE_PATTERN.search(str(units)):
            return str(name)
    for name in ds.variables:
        if str(name) in _TIME_NAMES:
            return str(name)
    for name, var in ds.variables.items():
        if var.attrs.get("axis") == "T":
            return str(name)
    return None


def configure_keys_for_input(options: Options) -> None:
    """Set y_key from the input_choice dispatch and auto-detect lat/lon/time.

    Opens the first matching file in the dataset's data folder and inspects
    its variable attributes (CF-style cascade) to identify the coordinate
    names. CLI overrides via --lat-var/--lon-var/--time-var take precedence.
    Time falls back to "time" when no time coord exists (e.g. AQUA_MODIS L3m,
    where the loader synthesises a "time" dim in its preprocess step).

    Args:
        options: Options object. Mutates time_key, y_key, lat_key, lon_key.

    Raises:
        KeyError: If options.input_choice is not in _INPUT_DISPATCH.
        FileNotFoundError: If the data folder contains no matching files.
        ValueError: If lat or lon cannot be detected and no override is given.
    """
    y_key, folder_attr, glob_pattern = _INPUT_DISPATCH[options.input_choice]
    options.y_key = y_key

    folder: Path = getattr(options, folder_attr)
    candidates = sorted(folder.glob(glob_pattern))
    if not candidates:
        raise FileNotFoundError(f"No files matching {glob_pattern!r} found in {folder}")
    representative = candidates[0]

    args = options.args
    with xr.open_dataset(representative) as ds:
        lat = args.lat_var or _detect_latitude(ds)
        lon = args.lon_var or _detect_longitude(ds)
        time = args.time_var or _detect_time(ds) or "time"

    if lat is None:
        raise ValueError(
            f"Could not detect latitude variable in {representative}; "
            f"specify it explicitly with --lat-var."
        )
    if lon is None:
        raise ValueError(
            f"Could not detect longitude variable in {representative}; "
            f"specify it explicitly with --lon-var."
        )

    options.lat_key = lat
    options.lon_key = lon
    options.time_key = time

    # Capture the inputs derive_spatial_chunks needs while the file is open.
    # Tolerate bogus overrides — leave defaults so callers still see a value
    # even if the resolved lat/lon name doesn't exist in the file.
    with xr.open_dataset(representative) as ds:
        try:
            options.n_lat_raw = int(ds.sizes[lat])
            options.n_lon_raw = int(ds.sizes[lon])
            options.on_disk_chunks = _extract_lat_lon_disk_chunks(
                ds[y_key], lat_key=lat, lon_key=lon
            )
        except KeyError:
            options.n_lat_raw = 0
            options.n_lon_raw = 0
            options.on_disk_chunks = None

    logging.info(
        f"Detected coords for {options.input_choice}: "
        f"time={time} lat={lat} lon={lon}, y={y_key}"
    )
    logging.info(
        f"First-file shape: {options.n_lat_raw} × {options.n_lon_raw}, "
        f"on-disk chunks: {options.on_disk_chunks}"
    )


def _extract_lat_lon_disk_chunks(
    da: xr.DataArray, *, lat_key: str, lon_key: str
) -> tuple[int, int] | None:
    """Pull the (lat, lon) entries out of a DataArray's HDF5 chunksizes.

    Args:
        da: The data variable, with .encoding["chunksizes"] possibly set.
        lat_key: Name of the latitude dim.
        lon_key: Name of the longitude dim.

    Returns:
        (cy, cx) — chunk size along (lat, lon) — or None if the encoding
        doesn't expose chunksizes, the dims don't include both axes, or the
        chunk along that axis equals the full dim (no useful chunking).
    """
    chunksizes = da.encoding.get("chunksizes")
    if chunksizes is None:
        return None
    if lat_key not in da.dims or lon_key not in da.dims:
        return None
    try:
        cy = int(chunksizes[da.dims.index(lat_key)])
        cx = int(chunksizes[da.dims.index(lon_key)])
    except (IndexError, TypeError):
        return None
    if cy >= int(da.sizes[lat_key]) and cx >= int(da.sizes[lon_key]):
        return None
    return (cy, cx)


def _calculate_spatial_chunks(
    n_files: int,
    n_lat_raw: int,
    n_lon_raw: int,
    xskip: int,
    on_disk_chunks: tuple[int, int] | None,
    block_bytes: float,
    max_graph_chunks: int,
    itemsize: int = 4,
    manual_lat: int | None = None,
    manual_lon: int | None = None,
) -> tuple[int, int, int]:
    """Pure spatial-chunk derivation for ``xr.open_mfdataset``.

    Implements the algorithm in
    /home/claudeuser/.claude/plans/you-just-modified-timeseries2color-py-ethereal-sunset.md

    Args:
        n_files: Number of files in the multi-file dataset.
        n_lat_raw: Latitude dim of one file before xskip.
        n_lon_raw: Longitude dim of one file before xskip.
        xskip: Stride used in the loader's preprocess isel.
        on_disk_chunks: HDF5 chunk shape (cy, cx) along (lat, lon), or None.
        block_bytes: Target per-tile RAM (the same as the existing option).
        max_graph_chunks: Hard cap on total dask chunks across all files.
        itemsize: Bytes per element after preprocess (float32 → 4).
        manual_lat: User override; takes precedence when both are given.
        manual_lon: User override; takes precedence when both are given.

    Returns:
        (lat_chunk, lon_chunk, tile_side).  When n_files <= 1 the chunks are
        -1 (let xarray fall back to one chunk per file, no graph hazard).
    """
    # Step 0 — single-file short-circuit.
    if n_files <= 1:
        return -1, -1, 1

    # Step 1 — post-xskip dims.
    n_lat = -(-n_lat_raw // xskip)
    n_lon = -(-n_lon_raw // xskip)

    # Step 2 — pre-load tile_side estimate.
    tile_side = max(
        1,
        min(
            n_lat,
            n_lon,
            int(np.sqrt(block_bytes / (n_files * itemsize))),
        ),
    )

    # Step 8 — manual override wins over derivation.
    if manual_lat is not None and manual_lon is not None:
        return (
            max(1, min(manual_lat, n_lat)),
            max(1, min(manual_lon, n_lon)),
            tile_side,
        )

    # Step 3 — graph-budget lower bound. Start at the continuous
    # approximation and bump up until ceiling effects no longer push the
    # exact chunk count over the budget.
    c_budget = int(np.ceil(np.sqrt(n_files * n_lat * n_lon / max_graph_chunks)))
    cap = min(n_lat, n_lon)
    while (
        c_budget < cap
        and n_files * -(-n_lat // c_budget) * -(-n_lon // c_budget)
        > max_graph_chunks
    ):
        c_budget += 1

    # Step 4 — combine bounds.
    c_lo = max(c_budget, tile_side)

    # Step 5 — apply per-axis cap.
    lat_chunk = min(c_lo, n_lat)
    lon_chunk = min(c_lo, n_lon)

    # Step 6 — disk alignment (only when xskip == 1 and on-disk chunks known).
    if xskip == 1 and on_disk_chunks is not None:
        cy, cx = on_disk_chunks
        lat_chunk = min(n_lat, cy * -(-lat_chunk // cy))
        lon_chunk = min(n_lon, cx * -(-lon_chunk // cx))

    # Step 7 — exact verify, fall back to one chunk per file if over budget.
    # (After Step 3's iteration this only fires when n_files alone > B.)
    n_total = n_files * -(-n_lat // lat_chunk) * -(-n_lon // lon_chunk)
    if n_total > max_graph_chunks:
        lat_chunk = n_lat
        lon_chunk = n_lon

    return lat_chunk, lon_chunk, tile_side


def derive_spatial_chunks(options: "Options") -> None:
    """Set ``options.lat_chunk`` / ``lon_chunk`` / ``tile_side`` from the
    actual dataset. Logs the rationale at INFO level.

    Must run after ``configure_keys_for_input`` (which populates
    ``n_lat_raw``, ``n_lon_raw``, ``on_disk_chunks``) and after the CLI
    backfill (which sets the final ``xskip``), but before
    ``load_input_data``.

    Args:
        options: Options object. Mutates ``lat_chunk``, ``lon_chunk``,
                 ``tile_side``.
    """
    # Recover n_files from the dispatch table without re-opening any file.
    _y_key, folder_attr, glob_pattern = _INPUT_DISPATCH[options.input_choice]
    folder: Path = getattr(options, folder_attr)
    candidates = sorted(folder.glob(glob_pattern))
    n_files = len(candidates[:: options.tskip])

    manual_lat = getattr(options.args, "lat_chunk", None)
    manual_lon = getattr(options.args, "lon_chunk", None)

    lat_chunk, lon_chunk, tile_side = _calculate_spatial_chunks(
        n_files=n_files,
        n_lat_raw=options.n_lat_raw,
        n_lon_raw=options.n_lon_raw,
        xskip=options.xskip,
        on_disk_chunks=options.on_disk_chunks,
        block_bytes=options.block_bytes,
        max_graph_chunks=options.max_graph_chunks,
        itemsize=4,
        manual_lat=manual_lat,
        manual_lon=manual_lon,
    )

    options.lat_chunk = lat_chunk
    options.lon_chunk = lon_chunk
    options.tile_side = tile_side

    n_lat_post = -(-options.n_lat_raw // options.xskip)
    n_lon_post = -(-options.n_lon_raw // options.xskip)
    if lat_chunk > 0 and lon_chunk > 0:
        n_total = (
            n_files
            * -(-n_lat_post // lat_chunk)
            * -(-n_lon_post // lon_chunk)
        )
    else:
        n_total = n_files  # single-file or sentinel
    logging.info(
        f"CHUNK INPUTS: n_files={n_files}, "
        f"n_lat_raw={options.n_lat_raw}, n_lon_raw={options.n_lon_raw}, "
        f"xskip={options.xskip}, on_disk={options.on_disk_chunks}"
    )
    logging.info(
        f"CHUNK FINAL: lat_chunk={lat_chunk}, lon_chunk={lon_chunk}, "
        f"tile_side={tile_side}; total dask chunks ≈ {n_total:,} "
        f"(budget {options.max_graph_chunks:,})"
    )
    if n_total > options.max_graph_chunks:
        logging.warning(
            f"CHUNK BUDGET: n_files={n_files} alone exceeds "
            f"max_graph_chunks={options.max_graph_chunks}; consider raising "
            f"--max-graph-chunks or increasing --tskip."
        )


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def parse_arguments(options: Options) -> None:
    """Parse command-line arguments.

    Args:
        options: Options object to store parsed arguments. Contains:
                     - my_name:  Name of the program.
                     - log_mode: Logging mode (default is logging.INFO).
                     - args:     Parsed arguments will be stored here.

    Returns:
        None, but updates options.args with parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=f"Spectral color mapping. {options.my_name} version {__version__}"
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
    )
    parser.add_argument(
        "description",
        type=str,
        nargs="?",
        default="",
        help="A description string encapsulated in quotes.",
    )
    parser.add_argument(
        "--input-choice",
        type=str,
        default=None,
        choices=options.input_choices_list,
        help=f"Data source (default: {options.input_choice}).",
    )
    parser.add_argument(
        "--lat-var",
        type=str,
        default=None,
        help="Override auto-detected latitude variable name.",
    )
    parser.add_argument(
        "--lon-var",
        type=str,
        default=None,
        help="Override auto-detected longitude variable name.",
    )
    parser.add_argument(
        "--time-var",
        type=str,
        default=None,
        help="Override auto-detected time variable name.",
    )
    parser.add_argument(
        "--xskip",
        type=int,
        default=None,
        help="Skip every N points in lat/lon (default: per-dataset).",
    )
    parser.add_argument(
        "--min-period",
        type=float,
        default=None,
        help="Minimum period in days (default: per-dataset).",
    )
    parser.add_argument(
        "--max-period",
        type=float,
        default=None,
        help="Maximum period in days (default: per-dataset).",
    )
    parser.add_argument(
        "--chunk-bytes",
        type=float,
        default=None,
        help=f"Target bytes per dask time chunk (default: {options.chunk_bytes:.0e}).",
    )
    parser.add_argument(
        "--block-bytes",
        type=float,
        default=None,
        help=f"Target bytes per loop lat/lon block (default: {options.block_bytes:.0e}).",
    )
    parser.add_argument(
        "--max-graph-chunks",
        type=int,
        default=None,
        help=(
            "Hard cap on total dask chunks across all opened files "
            f"(default: {options.max_graph_chunks:,}). Lower if open_mfdataset "
            "OOMs; raise for finer-grained spatial chunks."
        ),
    )
    parser.add_argument(
        "--lat-chunk",
        type=int,
        default=None,
        help=(
            "Override auto-derived lat chunk size for open_mfdataset. "
            "Must be paired with --lon-chunk."
        ),
    )
    parser.add_argument(
        "--lon-chunk",
        type=int,
        default=None,
        help=(
            "Override auto-derived lon chunk size for open_mfdataset. "
            "Must be paired with --lat-chunk."
        ),
    )
    parser.add_argument(
        "--use-kerchunk",
        action="store_true",
        help=(
            "Load multi-file datasets through a kerchunk virtual-zarr "
            "reference. Faster for large xskip values because zarr's "
            "chunk-aware indexing skips disk chunks that contain no "
            "target indices. Builds a sidecar .kerchunk.json on first run."
        ),
    )
    parser.add_argument(
        "--kerchunk-path",
        type=Path,
        default=None,
        help=(
            "Override the kerchunk reference file path "
            "(default: <data folder>/.kerchunk.json)."
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=options.dpi,
        help=f"DPI for output images (default: {options.dpi}).",
    )
    parser.add_argument(
        "-d", "--debug", action="store_true", help="Enable debug logging."
    )
    parser.add_argument(
        "--no-gmt", action="store_true", help="Skip GMT plot generation."
    )
    parser.add_argument(
        "--projection",
        type=int,
        choices=[1, 2, 3, 4],
        default=None,
        help="Map projection: 1=Robinson, 2=Winkel Tripel, 3=Mollweide, 4=Miller.",
    )
    parser.add_argument(
        "--use-old-funcs",
        action="store_true",
        help="Use Hughes & Williams 2010 functions instead of colour-science.",
    )
    parser.add_argument(
        "--land-color",
        type=str,
        default=None,
        help="Land color for matplotlib maps (default: '#333333').",
    )
    parser.add_argument(
        "--ocean-color",
        type=str,
        default=None,
        help="Ocean color for matplotlib maps (default: 'black').",
    )
    parser.add_argument(
        "--light-mode",
        action="store_true",
        help="Use light/white background instead of dark.",
    )
    parser.add_argument(
        "--borders", action="store_true", help="Show political boundaries on maps."
    )
    parser.add_argument(
        "--rivers", action="store_true", help="Show major rivers on maps."
    )
    parser.add_argument(
        "--coastline-resolution",
        type=str,
        default=None,
        choices=["110m", "50m"],
        help="Coastline resolution (default: '110m').",
    )
    parser.add_argument(
        "--fit-terms",
        type=str,
        nargs="+",
        default=None,
        choices=["constant", "trend", "accel", "annual", "semiannual"],
        help="Detrending terms (default: constant trend accel annual semiannual).",
    )
    options.args = parser.parse_args()
    if options.args.debug:
        options.log_mode = logging.DEBUG
    if options.args.input_choice is not None:
        options.input_choice = options.args.input_choice
    if options.args.no_gmt:
        options.plot_options.no_gmt_plots = 1
    if options.args.projection is not None:
        options.plot_options.projection = options.args.projection
    if options.args.use_old_funcs:
        options.use_new_funcs = False
    if options.args.land_color is not None:
        options.plot_options.land_color = options.args.land_color
    if options.args.ocean_color is not None:
        options.plot_options.ocean_color = options.args.ocean_color
    if options.args.light_mode:
        options.plot_options.dark_mode = False
    if options.args.borders:
        options.plot_options.show_borders = True
    if options.args.rivers:
        options.plot_options.show_rivers = True
    if options.args.coastline_resolution is not None:
        options.plot_options.coastline_resolution = options.args.coastline_resolution
    if options.args.chunk_bytes is not None:
        options.chunk_bytes = options.args.chunk_bytes
    if options.args.block_bytes is not None:
        options.block_bytes = options.args.block_bytes
    if options.args.max_graph_chunks is not None:
        options.max_graph_chunks = options.args.max_graph_chunks
    if options.args.fit_terms is not None:
        options.fit_terms = options.args.fit_terms
    if options.args.use_kerchunk:
        options.use_kerchunk = True
    if options.args.kerchunk_path is not None:
        options.kerchunk_path = options.args.kerchunk_path
    options.dpi = options.args.dpi


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main() -> None:
    """Main function."""
    start_time = dt.datetime.now()

    options = Options()
    parse_arguments(options)

    # Logging setup (replaces deleted logging_setup function)
    logging.basicConfig(
        level=options.log_mode,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Create output folder
    date_str = start_time.strftime("%Y%m%d-%H%M%S")
    options.output_folder = options.output_base / date_str
    if options.args.description:
        options.output_folder = options.output_folder.parent / (
            options.output_folder.name + " - " + options.args.description
        )
    options.output_folder.mkdir(parents=True, exist_ok=True)
    if not options.output_folder.is_dir():
        raise ValueError(f"!!! Problem creating {os.fspath(options.output_folder)}")

    # Copy support files
    for file_name in options.files_to_copy:
        shutil.copy(file_name, options.output_folder)

    # Zip the current script into the output folder
    zip_script(options)

    # Matplotlib setup
    plt.style.use("dark_background" if options.plot_options.dark_mode else "default")
    plt.rcParams["font.size"] = 14
    plt.rcParams["axes.linewidth"] = 2

    funcs_desc = " NEW functions" if options.use_new_funcs else " OLD functions"

    # Configure data keys per input_choice
    configure_keys_for_input(options)

    # Backfill None CLI values with per-dataset defaults
    for attr in ("xskip", "min_period", "max_period"):
        if getattr(options.args, attr) is None:
            setattr(options.args, attr, getattr(options, attr))
        else:
            setattr(options, attr, getattr(options.args, attr))

    # Derive (lat_chunk, lon_chunk) from dataset properties before the loader
    # is called.  This sets options.lat_chunk and options.lon_chunk, which the
    # multi-file loaders pass to xr.open_mfdataset.
    derive_spatial_chunks(options)

    # Coarse initial chunk_size for the initial open_mfdataset call.
    # Actual value is derived from options.chunk_bytes after we know grid dimensions.
    options.chunk_size = 100
    logging.info(f"Loading {options.input_choice}")
    input_data = load_input_data(options)

    # Derive final chunk_size and tile_side from target byte budgets
    _da = input_data[options.y_key]
    _itemsize = _da.dtype.itemsize
    _n_lat = _da.sizes[options.lat_key]
    _n_lon = _da.sizes[options.lon_key]
    _n_time = _da.sizes[options.time_key]
    _n_points = _n_lat * _n_lon

    options.chunk_size = max(
        1,
        min(
            _n_time,
            int(options.chunk_bytes / (_n_lat * _n_lon * _itemsize)),
        ),
    )
    # Square-ish (lat, lon) tile so each per-block read intersects only a
    # bounded set of on-disk chunks.  Cap by the smaller spatial dim.
    options.tile_side = max(
        1,
        min(
            _n_lat,
            _n_lon,
            int(np.sqrt(options.block_bytes / (_n_time * _itemsize))),
        ),
    )

    logging.info(
        f"DERIVED chunk_size = {options.chunk_size} "
        f"(target {options.chunk_bytes / 1e6:.0f} MB per chunk, "
        f"actual ≈ {options.chunk_size * _n_lat * _n_lon * _itemsize / 1e6:.0f} MB)"
    )
    logging.info(
        f"DERIVED tile_side = {options.tile_side} "
        f"(target {options.block_bytes / 1e6:.0f} MB per tile, "
        f"actual ≈ {_n_time * options.tile_side * options.tile_side * _itemsize / 1e6:.0f} MB)"
    )
    del _da

    # Rechunk along time only.  Spatial rechunking would produce a
    # dask graph with one task per (time_chunk × lat_tile × lon_tile),
    # which on a multi-thousand-file archive at 4320×8640 would itself
    # consume gigabytes of graph state.  Keeping lat/lon as one chunk
    # per timestep means each tile.values still streams chunk-by-chunk
    # under the synchronous scheduler, with peak working set ≈ one chunk.
    input_data = input_data.chunk({options.time_key: options.chunk_size})

    # --- DIAGNOSTIC: Tuning chunk_size and tile_side ---
    _da = input_data[options.y_key]
    _itemsize = _da.dtype.itemsize  # 4 for float32
    _n_time = _da.sizes[options.time_key]
    _n_lat = _da.sizes[options.lat_key]
    _n_lon = _da.sizes[options.lon_key]
    _n_points = _n_lat * _n_lon

    # Actual chunk sizes along time (may differ from requested if files contain 1 step each
    # and xarray can't fuse across files cleanly)
    _time_chunks = _da.chunks[_da.dims.index(options.time_key)]
    _effective_chunk = max(_time_chunks)  # largest per-chunk memory
    _n_chunks = len(_time_chunks)

    # Per-chunk memory footprint when dask materializes one chunk
    _chunk_bytes = _effective_chunk * _n_lat * _n_lon * _itemsize
    logging.info(
        f"TUNING chunk_size: requested={options.chunk_size}, "
        f"effective max along {options.time_key}={_effective_chunk} "
        f"({_n_chunks} chunks total)"
    )
    logging.info(
        f"TUNING chunk_size: peak per-chunk RAM ≈ "
        f"{_effective_chunk} × {_n_lat} × {_n_lon} × {_itemsize} B = "
        f"{_chunk_bytes / 1e9:.3f} GB"
    )

    if _chunk_bytes > 500e6:
        logging.warning(
            f"TUNING chunk_size: {_chunk_bytes / 1e9:.2f} GB per chunk is large "
        )
    elif _chunk_bytes < 10e6 and _n_chunks > 50:
        logging.warning(
            f"TUNING chunk_size: {_chunk_bytes / 1e6:.1f} MB per chunk is small and "
            f"there are {_n_chunks} chunks"
        )

    # Per-tile memory footprint in the main loop
    _tile_pts = options.tile_side * options.tile_side
    _tile_bytes = _n_time * _tile_pts * _itemsize
    _n_lat_tiles = (_n_lat + options.tile_side - 1) // options.tile_side
    _n_lon_tiles = (_n_lon + options.tile_side - 1) // options.tile_side
    _n_tiles = _n_lat_tiles * _n_lon_tiles
    logging.info(
        f"TUNING tile_side: {_n_tiles} tiles "
        f"({_n_lat_tiles} × {_n_lon_tiles}) of up to "
        f"{options.tile_side} × {options.tile_side} pixels each"
    )
    logging.info(
        f"TUNING tile_side: peak per-tile RAM ≈ "
        f"{_n_time} × {options.tile_side}² × {_itemsize} B = "
        f"{_tile_bytes / 1e9:.3f} GB"
    )
    if _tile_bytes > 1e9:
        logging.warning(
            f"TUNING tile_side: {_tile_bytes / 1e9:.2f} GB per tile is large "
        )
    elif _tile_bytes < 20e6 and _n_tiles > 20:
        logging.warning(
            f"TUNING tile_side: {_tile_bytes / 1e6:.1f} MB per tile is small and "
            f"there are {_n_tiles} tiles"
        )

    # Combined peak estimate (tile being built + one chunk live in dask's hands)
    _peak = _tile_bytes + _chunk_bytes
    logging.info(
        f"TUNING combined: rough peak RAM during loop ≈ "
        f"{_peak / 1e9:.3f} GB (tile + one chunk)"
    )
    del _da

    # --- DIAGNOSTIC: Check loaded data ---
    if 0:
        logging.info(f"DIAGNOSTIC: chunks = {input_data[options.y_key].chunks}")
        logging.info(f"DIAGNOSTIC: input_data dims = {dict(input_data.sizes)}")
        logging.info(f"DIAGNOSTIC: time dtype = {input_data[options.time_key].dtype}")
        _data_var = input_data[options.y_key]
        _total = int(_data_var.size)
        _stats = xr.Dataset(
            {
                "nan": _data_var.isnull().sum(),
                "vmin": _data_var.min(skipna=True),
                "vmax": _data_var.max(skipna=True),
            }
        ).compute()
        _nan_count = int(_stats["nan"])
        _vmin = float(_stats["vmin"])
        _vmax = float(_stats["vmax"])
        logging.info(
            f"DIAGNOSTIC: {options.y_key} NaN={_nan_count}, valid={_total - _nan_count}, total={_total}"
        )
        logging.info(
            f"DIAGNOSTIC: {options.y_key} finite range = [{_vmin:.6e}, {_vmax:.6e}]"
        )
        if _vmax > 1e10:
            logging.error(
                f"DIAGNOSTIC: WARNING - {options.y_key} has extreme values (max={_vmax:.2e}), fill values may not be masked!"
            )

    if options.xskip > 1:
        logging.info(f"Selecting one lat/lon point in every {options.xskip**2} points.")
    else:
        logging.info("Using all lat/lon points.")

    # Ensure even number of time steps for NFFT
    if input_data.sizes[options.time_key] % 2 == 1:
        logging.info(
            f"Deleting last data point because number of time stamps needs to be even for NFFT. Before deletion: {input_data.sizes[options.time_key] = }"
        )
        input_data = input_data.isel({options.time_key: slice(None, -1)})

    # Validate min/max period against available NFFT periods
    # Pick a grid point with valid data (index 0 may be land/pole)
    _first_time = input_data[options.y_key].isel({options.time_key: 0}).values
    _valid_pts = np.argwhere(np.isfinite(_first_time))
    if len(_valid_pts) > 0:
        _mid = len(_valid_pts) // 2
        _vlat, _vlon = int(_valid_pts[_mid, 0]), int(_valid_pts[_mid, 1])
    else:
        _vlat, _vlon = 0, 0
    sliced_data = input_data.isel({options.lat_key: _vlat, options.lon_key: _vlon})
    power_spectrum = nfft_power(options, sliced_data)
    power_spectrum = convert_spectrum_from_frequency_to_period(power_spectrum)

    logging.info(
        f"{options.min_period = } and {np.min(power_spectrum.period.values) = }"
    )
    if options.min_period < np.min(
        power_spectrum.period.values
    ) or options.min_period >= np.max(power_spectrum.period.values):
        logging.error(
            f"!!! WARNING!!! originally {options.min_period = } but {np.min(power_spectrum.period.values) = } and {np.max(power_spectrum.period.values) = }"
        )
        options.min_period = np.min(power_spectrum.period.values)
        logging.error(
            f"So now {options.min_period = } which equals {np.min(power_spectrum.period.values) = }"
        )
    logging.info(
        f"{options.max_period = } and {np.max(power_spectrum.period.values) = }"
    )
    if options.max_period <= np.min(
        power_spectrum.period.values
    ) or options.max_period > np.max(power_spectrum.period.values):
        logging.error(
            f"!!! WARNING!!! originally {options.max_period = } but {np.min(power_spectrum.period.values) = } and {np.max(power_spectrum.period.values) = }"
        )
        options.max_period = np.max(power_spectrum.period.values)
        logging.error(
            f"So now {options.max_period = } which equals {np.max(power_spectrum.period.values) = }"
        )

    # Load CIE color matching functions
    cie = load_cie_functions(options)

    # ===================================================================
    # HOT PATH — optimized computation loop
    # ===================================================================
    hot_start = timeit.default_timer()

    # === PRE-COMPUTATION (outside loop) ===
    # No .stack here: stacking lat × lon eagerly builds a pandas.MultiIndex
    # of size n_lat × n_lon (~37M entries on AQUA_MODIS 4km), which costs
    # multiple GB of RAM by itself.  Iterate tiles in (lat, lon) directly.
    times = input_data[options.time_key].values
    lat_coord = input_data[options.lat_key]
    lon_coord = input_data[options.lon_key]
    data_lazy = input_data[options.y_key]  # dask-backed DataArray
    n_times = len(times)
    n_lat = data_lazy.sizes[options.lat_key]
    n_lon = data_lazy.sizes[options.lon_key]
    n_points = n_lat * n_lon

    logging.info(
        f"Full data_2d would be {n_times} × {n_points} × 4 bytes = "
        f"{n_times * n_points * 4 / 1e9:.2f} GB — loading in tiles."
    )

    # Time axis in days (shared across all grid points)
    time_days = (times - times[0]) / np.timedelta64(1, "D")

    # --- DIAGNOSTIC: Check time conversion ---
    logging.info(f"DIAGNOSTIC: times dtype = {times.dtype}")
    logging.info(
        f"DIAGNOSTIC: time_days range = [{time_days.min():.4f}, {time_days.max():.4f}] days"
    )
    if time_days.max() < 1.0:
        logging.error(
            "DIAGNOSTIC: CRITICAL - time_days span < 1 day! Time conversion may be wrong."
        )

    # Ensure even length for NFFT
    N = len(time_days)
    if N % 2:
        logging.info(f"Trimming last time step for even NFFT length: {N} -> {N - 1}")
        time_days = time_days[:-1]
        n_times -= 1
        N -= 1

    # Detrending design matrix (shared, controlled by options.fit_terms)
    design_matrix = build_design_matrix(time_days, options.fit_terms)

    # NFFT parameters (shared)
    x_min = time_days.min()
    x_range = time_days.max() - x_min
    x_norm = (time_days - x_min) / x_range - 0.5
    k = -(N // 2) + np.arange(N)
    xf = k / x_range
    xf_half = xf[N // 2 + 1 :]  # positive frequencies, ascending

    # Period array (ascending) for mapping power spectrum to wavelength
    periods_ascending = (1.0 / xf_half)[::-1]
    # PSD-to-SPD conversion factor: sigma^2 (Hughes & Williams 2010, eq. A11)
    freq_sq_ascending = (1.0 / periods_ascending) ** 2
    cie_wavelengths = cie["wavelength"].values
    n_wl = len(cie_wavelengths)
    wavelength_targets = np.linspace(options.min_period, options.max_period, n_wl)

    # CIE color matching functions for sd_to_XYZ
    cmfs = MSDS_CMFS_STANDARD_OBSERVER["CIE 1931 2 Degree Standard Observer"]

    # Pre-compute for old (Riemann sum) path
    cie_x_vals = cie["x"].values  # shape (n_wl,)
    cie_y_vals = cie["y"].values
    cie_z_vals = cie["z"].values
    wavelength_step = np.diff(cie_wavelengths).mean()

    A_old = A_OLD_XYZ_TO_SRGB

    # === MAIN LOOP — tile-wise (lat × lon) to keep memory bounded ===
    # Shape (n_lat, n_lon, 3) so per-pixel writes are local within a tile.
    # Reshaped to flat (n_points, 3) after the loop for the existing
    # post-processing.
    xyz_results = np.full((n_lat, n_lon, 3), np.nan)
    tile_side = options.tile_side

    with tqdm(total=n_points, desc="Grid points") as pbar:
        for lat0 in range(0, n_lat, tile_side):
            lat1 = min(lat0 + tile_side, n_lat)
            for lon0 in range(0, n_lon, tile_side):
                lon1 = min(lon0 + tile_side, n_lon)
                # Materialize just this tile (each .values triggers a fresh
                # dask compute that reads only the chunks intersecting the tile).
                tile = data_lazy.isel(
                    {
                        options.lat_key: slice(lat0, lat1),
                        options.lon_key: slice(lon0, lon1),
                    }
                ).values
                if N < tile.shape[0]:  # handle the even-length trim
                    tile = tile[:N, :, :]

                lat_t = tile.shape[1]
                lon_t = tile.shape[2]
                for jl in range(lat_t):
                    ilat = lat0 + jl
                    for jo in range(lon_t):
                        ilon = lon0 + jo
                        pbar.update(1)

                        y = tile[:, jl, jo]

                        # Skip all-NaN columns (land)
                        if np.all(np.isnan(y)):
                            continue

                        # Drop NaN samples — use only real data points
                        nan_idx = np.isnan(y)
                        if np.any(nan_idx):
                            valid_idx = ~nan_idx
                            n_valid_pts = int(valid_idx.sum())
                            # The nfft library's kernel window half-width m is
                            # estimated as 9 for the default tolerance, and
                            # the internal assertion m <= n // 2 (n = N * 3)
                            # fails when N < 6.
                            if n_valid_pts < 6:
                                continue
                            y = y[valid_idx]
                            y_time_days = time_days[valid_idx]
                            y_N = n_valid_pts
                            if y_N % 2:
                                y = y[:-1]
                                y_time_days = y_time_days[:-1]
                                y_N -= 1
                            y_x_min = y_time_days.min()
                            y_x_range = y_time_days.max() - y_x_min
                            y_x_norm = (y_time_days - y_x_min) / y_x_range - 0.5
                            y_design = build_design_matrix(
                                y_time_days, options.fit_terms
                            )
                            y_k = -(y_N // 2) + np.arange(y_N)
                            y_xf_half = (y_k / y_x_range)[y_N // 2 + 1 :]
                            y_periods_asc = (1.0 / y_xf_half)[::-1]
                            y_freq_sq_asc = (1.0 / y_periods_asc) ** 2
                        else:
                            y_N = N
                            y_x_norm = x_norm
                            y_design = design_matrix
                            y_periods_asc = periods_ascending
                            y_freq_sq_asc = freq_sq_ascending

                        # 1. Detrend with numpy lstsq
                        params, _, _, _ = np.linalg.lstsq(y_design, y, rcond=None)
                        detrended = y - y_design @ params

                        # 2. NFFT
                        f_k = nfft.nfft(y_x_norm, detrended)
                        power = np.abs(f_k) ** 2

                        # 3. Positive frequencies, reverse to ascending period
                        power_half = power[y_N // 2 + 1 :]
                        power_ascending = power_half[::-1]

                        # 3b. PSD -> SPD via sigma^2
                        # (Hughes & Williams 2010, eq. A11)
                        power_ascending = power_ascending * y_freq_sq_asc

                        # 4. Map power spectrum to wavelength grid
                        pt_periods_ext = y_periods_asc.copy()
                        if options.min_period not in y_periods_asc:
                            pt_periods_ext = np.append(
                                pt_periods_ext, options.min_period
                            )
                        if options.max_period not in y_periods_asc:
                            pt_periods_ext = np.append(
                                pt_periods_ext, options.max_period
                            )
                        pt_periods_ext = np.sort(pt_periods_ext)
                        pt_nearest_idx = np.array(
                            [
                                np.argmin(np.abs(y_periods_asc - p))
                                for p in pt_periods_ext
                            ]
                        )
                        pt_mask = (pt_periods_ext >= options.min_period) & (
                            pt_periods_ext <= options.max_period
                        )
                        pt_interp_periods = pt_periods_ext[pt_mask]
                        pt_interp_nearest = pt_nearest_idx[pt_mask]
                        power_interp_src = power_ascending[pt_interp_nearest]
                        mapped_power = np.interp(
                            wavelength_targets,
                            pt_interp_periods,
                            power_interp_src,
                        )

                        # 5. Spectrum to XYZ
                        if options.use_new_funcs:
                            spd = SpectralDistribution(mapped_power, cie_wavelengths)
                            xyz_results[ilat, ilon, :] = sd_to_XYZ(spd, cmfs)
                        else:
                            # Riemann sum (Hughes & Williams 2010)
                            xyz_results[ilat, ilon, 0] = (
                                mapped_power * cie_x_vals
                            ).sum() * wavelength_step
                            xyz_results[ilat, ilon, 1] = (
                                mapped_power * cie_y_vals
                            ).sum() * wavelength_step
                            xyz_results[ilat, ilon, 2] = (
                                mapped_power * cie_z_vals
                            ).sum() * wavelength_step

    # Cleanup — flatten for the existing post-loop code
    del input_data, data_lazy
    gc.collect()
    xyz_results = xyz_results.reshape(n_points, 3)

    logging.info("Main loop complete.")

    # --- DIAGNOSTIC: Check xyz_results ---
    _xyz_valid = int(np.any(np.isfinite(xyz_results), axis=1).sum())
    _xyz_nan = int(np.all(np.isnan(xyz_results), axis=1).sum())
    logging.info(
        f"DIAGNOSTIC: xyz_results: {_xyz_valid} valid, {_xyz_nan} all-NaN, {n_points} total"
    )
    if _xyz_valid > 0:
        _y_valid = xyz_results[np.isfinite(xyz_results[:, 1]), 1]
        logging.info(
            f"DIAGNOSTIC: Y range = [{_y_valid.min():.6e}, {_y_valid.max():.6e}]"
        )
    else:
        logging.error("DIAGNOSTIC: CRITICAL - ALL xyz_results are NaN!")

    # === POST-LOOP: Normalize XYZ ===
    max_y = np.nanmax(xyz_results[:, 1])
    logging.info(f"The highest value of 'y' is: {max_y}")
    xyz_results /= max_y

    # Raise Y to power (brightens dark areas)
    logging.info(
        f"Raising y to power {options.thepower} while keeping chromaticity constant. (This brightens dark areas.)"
    )
    factor = np.power(xyz_results[:, 1], 1.0 - options.thepower)
    xyz_results /= factor[:, np.newaxis]

    # === POST-LOOP: Vectorized XYZ -> sRGB (single call, replaces groupby.map) ===
    logging.info("Converting to RGB...")
    if options.use_new_funcs:
        RGB_results = XYZ_to_sRGB(xyz_results)  # shape (n_points, 3), includes gamma
    else:
        # Linear matrix multiply (Hughes & Williams 2010, eq. A2); gamma applied later
        RGB_results = (A_old @ xyz_results.T).T  # shape (n_points, 3)

    # === POST-LOOP: Vectorized fix_gamut (replaces groupby.map) ===
    logging.info("Fixing RGB out-of-gamut values and normalizing...")
    rgb_arr = RGB_results.T.copy()  # shape (3, n_points) for easier per-channel access
    del RGB_results  # Clear memory

    nan_mask = np.any(np.isnan(rgb_arr), axis=0)
    min_vals = np.full(n_points, 0.0)
    min_vals[~nan_mask] = rgb_arr[:, ~nan_mask].min(axis=0)
    needs_fix = (~nan_mask) & (min_vals < 0)

    if np.any(needs_fix):
        offset = -min_vals[needs_fix] + 1.0 / 255.0
        y_vals = xyz_results[needs_fix, 1]
        fix_factor = y_vals / (y_vals + offset)
        rgb_arr[:, needs_fix] = (rgb_arr[:, needs_fix] + offset) * fix_factor

    # Normalize RGB
    highest_value = np.nanmax(rgb_arr)
    rgb_arr /= highest_value

    # --- DIAGNOSTIC: Check normalization ---
    logging.info(
        f"DIAGNOSTIC: max_y = {max_y:.6e}, highest_value = {highest_value:.6e}"
    )
    _rgb_valid_mask = ~np.any(np.isnan(rgb_arr), axis=0)
    _n_rgb_valid = int(_rgb_valid_mask.sum())
    if _n_rgb_valid > 0:
        _rgb_valid_vals = rgb_arr[:, _rgb_valid_mask]
        logging.info(
            f"DIAGNOSTIC: RGB valid points: {_n_rgb_valid}, range=[{_rgb_valid_vals.min():.6f}, {_rgb_valid_vals.max():.6f}]"
        )
    else:
        logging.error("DIAGNOSTIC: CRITICAL - ALL RGB values are NaN!")

    # === POST-LOOP: Apply BT.709 gamma correction (Hughes & Williams 2010, eqs. A7-A8) ===
    if not options.use_new_funcs:
        gamma_inv = 0.45
        crit = 0.018
        h = 4.506813168
        g_coeff = -0.09914989
        f_coeff = 1.09914989
        for ch in range(3):
            vals = rgb_arr[ch]
            low = (~np.isnan(vals)) & (vals <= crit)
            high = (~np.isnan(vals)) & (vals > crit)
            vals[low] = vals[low] * h
            vals[high] = f_coeff * np.power(vals[high], gamma_inv) + g_coeff

    # === Assemble xarray Dataset directly on (lat, lon) — no unstack ===
    xyz_3d = xyz_results.reshape(n_lat, n_lon, 3)
    rgb_3d = rgb_arr.reshape(3, n_lat, n_lon)
    dims2d = (options.lat_key, options.lon_key)
    coords2d = {options.lat_key: lat_coord, options.lon_key: lon_coord}
    spectral_color_maps = xr.Dataset(
        {
            "x": (dims2d, xyz_3d[:, :, 0]),
            "y": (dims2d, xyz_3d[:, :, 1]),
            "z": (dims2d, xyz_3d[:, :, 2]),
            "red": (dims2d, rgb_3d[0]),
            "green": (dims2d, rgb_3d[1]),
            "blue": (dims2d, rgb_3d[2]),
        },
        coords=coords2d,
    )

    # Stop the timer
    hot_end = timeit.default_timer()
    time_taken = hot_end - hot_start
    logging.info(f"Time taken: {time_taken:.2f} seconds")

    # Populate dataclasses from actual computed data
    populate_results(options, spectral_color_maps)

    # ===================================================================
    # Save outputs and plot
    # ===================================================================

    map_projection = _get_cartopy_projection(options.plot_options.projection)

    # Plot composite RGB map
    plot_rgb_map(
        spectral_color_maps["red"],
        spectral_color_maps["green"],
        spectral_color_maps["blue"],
        lat_key=options.lat_key,
        lon_key=options.lon_key,
        projection=map_projection,
        title=f"{options.input_choice}, periods {options.min_period / 7.0:.0f}-{options.max_period / 7.0:.0f} weeks,{funcs_desc}",
        output_path=options.output_folder
        / f"{options.input_choice}_map_matplotlib.png",
        plot_options=options.plot_options,
    )

    if not options.plot_options.no_gmt_plots:
        # Save XYZ tristimulus NetCDFs
        for thekey in ["x", "y", "z"]:
            da = spectral_color_maps[thekey]
            filename = options.output_folder / f"{thekey}.nc"
            da.coords[options.lat_key].attrs = {
                "units": "degrees_north",
                "long_name": "Latitude",
            }
            da.coords[options.lon_key].attrs = {
                "units": "degrees_east",
                "long_name": "Longitude",
            }
            da.attrs["long_name"] = f"Spectral color {thekey}"
            da.attrs["units"] = "dimensionless"
            logging.info(f"Saving {os.fspath(filename)}...")
            da.to_netcdf(filename)

        # Save RGB NetCDFs at [0, 1] range (matches reference output)
        for thekey in ["red", "green", "blue"]:
            da = spectral_color_maps[thekey]
            filename = options.output_folder / f"{thekey}.nc"
            da.coords[options.lat_key].attrs = {
                "units": "degrees_north",
                "long_name": "Latitude",
            }
            da.coords[options.lon_key].attrs = {
                "units": "degrees_east",
                "long_name": "Longitude",
            }
            da.attrs["long_name"] = f"Spectral color {thekey}"
            da.attrs["units"] = "dimensionless"
            logging.info(f"Saving {os.fspath(filename)}...")
            da.to_netcdf(filename)

        # Save RGB NetCDFs as uint8 [0, 255] for GMT.  GMT's grdimage
        # three-grid RGB mode normalizes float grids independently per
        # channel, destroying color relationships.  Uint8 grids are
        # used as-is.  All files use variable name "z" (GMT's default)
        # and NETCDF4_CLASSIC format for maximum compatibility.
        for thekey in ["red", "green", "blue"]:
            uint8_values = _float01_to_gmt_uint8(spectral_color_maps[thekey].values)
            gmt_da = xr.DataArray(
                uint8_values,
                coords=spectral_color_maps[thekey].coords,
                dims=spectral_color_maps[thekey].dims,
                name="z",
            )
            gmt_file = options.output_folder / f"{thekey}_gmt.nc"
            gmt_da.coords[options.lat_key].attrs = {
                "units": "degrees_north",
                "long_name": "Latitude",
            }
            gmt_da.coords[options.lon_key].attrs = {
                "units": "degrees_east",
                "long_name": "Longitude",
            }
            gmt_da.attrs["long_name"] = f"Spectral color {thekey} (scaled 0-255)"
            gmt_da.attrs["units"] = "dimensionless"
            logging.info(f"Saving {os.fspath(gmt_file)} (uint8 for GMT)...")
            gmt_da.to_netcdf(
                gmt_file,
                encoding={"z": {"dtype": "uint8", "_FillValue": None}},
            )
        logging.info("Finished saving NetCDFs.")

        write_rgb_colorscale(options, cie)
        write_gmt_scripts(options)
        run_gmt_scripts(options)

    logging.info("All finished!")
    logging.info(
        "Download and analyze chlorophyll data, as well as sea surface salinity data!"
    )

    end_time = dt.datetime.now()
    elapsed_time = end_time - start_time
    logging.info(f"Elapsed time: {elapsed_time}")
    logging.shutdown()


# ---------------------------------------------------------------------------
# Pure helper functions (no options needed)
# ---------------------------------------------------------------------------


def build_design_matrix(time_days: np.ndarray, fit_terms: list[str]) -> np.ndarray:
    """Build an OLS design matrix from a list of fit terms.

    Args:
        time_days: Time axis in days (relative to start).
        fit_terms: List of terms from {"constant", "trend", "accel",
                   "annual", "semiannual"}.  "annual" and "semiannual"
                   each contribute a sin/cos pair (2 columns).

    Returns:
        2D array of shape (len(time_days), n_columns).
    """
    columns: list[np.ndarray] = []
    for term in fit_terms:
        if term == "constant":
            columns.append(np.ones_like(time_days))
        elif term == "trend":
            columns.append(time_days)
        elif term == "accel":
            columns.append(time_days**2)
        elif term == "annual":
            columns.append(np.sin(2 * np.pi * time_days / DAYS_IN_YEAR))
            columns.append(np.cos(2 * np.pi * time_days / DAYS_IN_YEAR))
        elif term == "semiannual":
            columns.append(np.sin(4 * np.pi * time_days / DAYS_IN_YEAR))
            columns.append(np.cos(4 * np.pi * time_days / DAYS_IN_YEAR))
    return np.column_stack(columns)


def _float01_to_gmt_uint8(arr: np.ndarray) -> np.ndarray:
    """Convert a float [0, 1] array to uint8 [0, 255] for GMT grdimage.

    GMT's three-grid RGB mode normalizes float grids independently per channel,
    destroying color relationships.  Uint8 grids are used as-is.

    NaN values are replaced with 0 (black) because uint8 cannot represent NaN.
    Values are clipped to [0, 255] before casting to handle floating-point
    overshoot from gamma correction.

    Args:
        arr: 2D numpy array with values in [0, 1], possibly containing NaN.

    Returns:
        2D numpy array of dtype uint8 with values in [0, 255].
    """
    scaled = np.nan_to_num(arr, nan=0.0) * 255.0
    return np.clip(scaled, 0.0, 255.0).astype(np.uint8)


# ---------------------------------------------------------------------------
# Populate dataclasses from computed results
# ---------------------------------------------------------------------------


def populate_results(options: Options, spectral_color_maps: xr.Dataset) -> None:
    """Populate options.results and options.plot_options from actual computed data.

    Called after the hot path completes and spectral_color_maps is assembled.
    Replaces the old init_results/init_plot_options dummy values with real
    lat/lon coordinates, titles, units, and spatial extent.

    Args:
        options:              Options object (results, plot_options, grid mutated).
        spectral_color_maps:  xarray Dataset with x, y, z, red, green, blue
                              variables on (lat_key, lon_key) coordinates.
    """
    results = options.results
    lats = spectral_color_maps.coords[options.lat_key].values
    lons = spectral_color_maps.coords[options.lon_key].values

    # Populate latlon from actual computed grid
    results.latlon.lat = lats
    results.latlon.lon = lons
    results.latlon.outputs = [
        spectral_color_maps["red"].values,
        spectral_color_maps["green"].values,
        spectral_color_maps["blue"].values,
    ]

    # Spatial extent from actual data
    results.minlat = float(lats.min())
    results.maxlat = float(lats.max())
    results.minlon = float(lons.min())
    results.maxlon = float(lons.max())

    # Real metadata
    results.titles = ["Spectral Colors"]
    results.units = ["days"]
    results.rgb = [1.0, 1.0, 1.0]  # 3 elements signals RGB mode

    # Populate grid
    options.grid.lat = lats.tolist()
    options.grid.lon = lons.tolist()
    options.grid.maxlat = results.maxlat
    options.grid.minlat = results.minlat
    options.grid.maxlon = results.maxlon
    options.grid.minlon = results.minlon

    # Plot options
    options.plot_options.outputfolder = os.fspath(options.output_folder) + os.sep
    options.plot_options.just_the_filenames = [
        "red_gmt.nc",
        "green_gmt.nc",
        "blue_gmt.nc",
    ]
    options.plot_options.dpi = options.dpi


# ---------------------------------------------------------------------------
# GMT pure helpers (no options needed)
# ---------------------------------------------------------------------------


def write_gmt_defs(new_fp: BinaryIO) -> None:
    """Write clarifying GMT definitions to the script file."""
    new_fp.write(b"#######################################################\n")
    new_fp.write(b"#Clarifying definitions. Do not change!################\n")
    new_fp.write(b'start=" -K " #Should always redirect using > to write new PS.\n')
    new_fp.write(
        b'middle=" -O -K " #Should always redirect using >> to append to PS.\n'
    )
    new_fp.write(b'end=" -O " #Should always redirect using >> to append to PS.\n')
    new_fp.write(b"#######################################################\n")


def write_gmt_coastlines(new_fp: BinaryIO) -> None:
    """Write GMT coastlines commands to the script file."""
    new_fp.write(
        b"gmt pscoast -W$coast_thk/$coast_color $coast_res $range $projection $map_pos $middle >> $plot_base.ps\n"
    )
    new_fp.write(b"if [ $coastlines == 2 ]\n")
    new_fp.write(b"then\n")
    new_fp.write(
        b"  gmt psxy -N $coast_file -: -Sc$coast_thk -W$coast_thk/$coast_color $range $projection $map_pos $middle >> $plot_base.ps\n"
    )
    new_fp.write(b"fi\n")


def finish_flip_backgrounds(flip_fp: BinaryIO) -> None:
    """Finish writing the flip_backgrounds.sh script."""
    flip_fp.write(b"for current_base in $all_bases\n")
    flip_fp.write(b"do\n")
    flip_fp.write(b"  #Change black to cyan temporarily.\n")
    flip_fp.write(
        b"  convert $current_base.png -fill cyan -opaque black $current_base.png\n"
    )
    flip_fp.write(b"  #Change white to black.\n")
    flip_fp.write(
        b"  convert $current_base.png -fill black -opaque white $current_base.png\n"
    )
    flip_fp.write(b"  #Change temporary cyan to white.\n")
    flip_fp.write(
        b"  convert $current_base.png -fill white -opaque cyan $current_base.png\n"
    )
    flip_fp.write(b"done\n")


def finish_trim(trim_fp: BinaryIO) -> None:
    """Finish writing the trim.sh script."""
    trim_fp.write(b"for current_base in $all_bases\n")
    trim_fp.write(b"do\n")
    trim_fp.write(b"  convert $current_base.png -trim $current_base.png\n")
    trim_fp.write(b"done\n")


# ---------------------------------------------------------------------------
# Data / science functions
# ---------------------------------------------------------------------------


def load_cie_functions(options: Options) -> xr.Dataset:
    """Load CIE 1931 color matching functions from CSV.

    Args:
        options: Options object. Contains:
                     - cie_file: Path to the CIE CSV file.

    Returns:
        Dataset with 'x', 'y', 'z' variables and 'wavelength' coordinate.
    """
    cie_path = Path(options.cie_file)

    data: dict[str, list[float]] = {"x": [], "y": [], "z": []}
    wavelengths: list[float] = []

    try:
        with cie_path.open("r") as in_fp:
            reader = csv.reader(in_fp)
            for row in reader:
                wavelengths.append(float(row[0]))
                for i, val in enumerate(row[1:], start=0):
                    data[list(data.keys())[i]].append(float(val))
    except OSError:
        logging.error(f"The CIE file, {os.fspath(cie_path)}, failed to open.")

    da = xr.Dataset(
        {var: ("wavelength", data[var]) for var in data},
        coords={"wavelength": wavelengths},
    )
    da.attrs["title"] = "CIE 1931 color matching functions"
    da.coords["wavelength"].attrs["units"] = "nm"
    return da


def nfft_power(options: Options, timeseries: xr.Dataset) -> xr.Dataset:
    """Compute the power spectrum of a timeseries using NFFT.

    Args:
        options:    Options object. Contains:
                        - time_key: Time coordinate name.
                        - y_key: Variable name.
        timeseries: Input timeseries dataset.

    Returns:
        Dataset with 'power' variable and 'frequency' coordinate.
    """
    x = (
        timeseries[options.time_key] - timeseries[options.time_key][0]
    ).values.squeeze() / np.timedelta64(1, "D")
    y = timeseries[options.y_key].values.squeeze()

    # If timeseries has an odd number of points,
    # remove the last data point, then calculate min, range.
    N = 1
    while N % 2:
        N = len(x)
        x_min = np.min(x)
        x_range = np.max(x) - np.min(x)
        x_norm = (x - x_min) / x_range - 0.5
        if N % 2:
            logging.error(
                f"!!! WARNING!!! LENGTH NEEDS TO BE EVEN FOR NFFT, BUT: {len(x) = }"
            )
            logging.error("!!! DELETING LAST DATA POINT!")
            x = np.delete(x, -1)
            y = np.delete(y, -1)

    # Define Fourier modes
    k = -(N // 2) + np.arange(N)
    # Convert Fourier modes to frequencies
    xf = k / x_range

    # Perform NFFT.
    f_k = nfft.nfft(x_norm, y)

    # Compute power spectrum, which is the square of the absolute value of the Fourier Transform
    power_spectrum = np.abs(f_k) ** 2

    # Only take the positive frequencies. Since the output is symmetric, this will not lose any information.
    power_spectrum = power_spectrum[N // 2 + 1 :]
    xf_half = xf[N // 2 + 1 :]

    # Create xarray DataArray with coordinates
    power_spectrum_da = xr.DataArray(
        power_spectrum, coords=[("frequency", xf_half)], name="power"
    )

    # Convert this DataArray to a Dataset
    spectrum = power_spectrum_da.to_dataset()
    spectrum["frequency"].attrs["units"] = "1/days"

    return spectrum


def convert_spectrum_from_frequency_to_period(spectrum: xr.Dataset) -> xr.Dataset:
    """Convert a spectrum's coordinate from frequency to period.

    Args:
        spectrum: Dataset with 'frequency' coordinate and 'power' variable.

    Returns:
        Dataset with 'period' coordinate, sorted by ascending period.
    """
    freq_units = spectrum["frequency"].attrs.get("units", None)

    units_map = {"1/days": "days", "Hz": "seconds", "1/years": "years"}

    period = 1.0 / spectrum["frequency"]

    if freq_units in units_map:
        period.attrs["units"] = units_map[freq_units]
    else:
        logging.error(f"!!!WARNING!!! DID NOT RECOGNIZE {freq_units = }")

    new_spectrum = spectrum.assign_coords(period=("frequency", period.data))
    new_spectrum = new_spectrum.swap_dims({"frequency": "period"}).drop_vars(
        "frequency"
    )
    new_spectrum = new_spectrum.sortby("period")

    return new_spectrum


# ---------------------------------------------------------------------------
# Plot functions
# ---------------------------------------------------------------------------

_CARTOPY_PROJECTIONS: dict[int, type[ccrs.Projection]] = {
    1: ccrs.Robinson,
    2: ccrs.Robinson,  # Winkel Tripel not available in cartopy
    3: ccrs.Mollweide,
    4: ccrs.Miller,
}


def _get_cartopy_projection(projection_code: int) -> ccrs.Projection:
    """Return a cartopy projection matching the given integer code.

    Args:
        projection_code: Projection choice (1=Robinson, 2=Winkel Tripel
                         fallback to Robinson, 3=Mollweide, 4=Miller).

    Returns:
        A cartopy CRS Projection instance with central_longitude=180
        (appropriate for 0-360 longitude data).
    """
    if projection_code == 2:
        logging.warning(
            "Winkel Tripel is not available in cartopy. Using Robinson as fallback."
        )
    proj_cls = _CARTOPY_PROJECTIONS.get(projection_code, ccrs.Robinson)
    return proj_cls(central_longitude=180)


def _add_map_features(
    ax: Any,  # noqa: ANN401 — cartopy GeoAxes; no type stubs available
    plot_options: PlotOptions,
) -> None:
    """Add geographic features (land, ocean, coastlines, etc.) to a cartopy axis.

    Args:
        ax:           A GeoAxes with a cartopy projection.
        plot_options: PlotOptions controlling which features to draw.
    """
    coastline_color = "darkgray" if plot_options.dark_mode else "black"
    border_color = "gray" if plot_options.dark_mode else "#666666"
    label_color = "white" if plot_options.dark_mode else "black"

    # ax.add_feature(cfeature.OCEAN, facecolor=plot_options.ocean_color, zorder=0)
    # ax.add_feature(
    #     cfeature.LAND,
    #     facecolor=plot_options.land_color,
    #     zorder=1,
    #     edgecolor="none",
    # )
    ax.coastlines(
        resolution=plot_options.coastline_resolution,
        color=coastline_color,
        linewidth=0.1,
        zorder=3,
    )
    if plot_options.show_borders:
        ax.add_feature(
            cfeature.BORDERS,
            linewidth=0.3,
            edgecolor=border_color,
            zorder=4,
        )
    if plot_options.show_rivers:
        ax.add_feature(
            cfeature.RIVERS,
            linewidth=0.3,
            edgecolor="steelblue",
            zorder=4,
        )

    # Gridline labels only work on certain projections (PlateCarree, Mercator,
    # Miller).  For others cartopy silently ignores draw_labels, so it is safe
    # to always request them.
    gl = ax.gridlines(
        draw_labels=True,
        color=coastline_color,
        alpha=0.25,
        linestyle="--",
        linewidth=0.4,
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"color": label_color}
    gl.ylabel_style = {"color": label_color}


def plot_rgb_map(
    red: xr.DataArray,
    green: xr.DataArray,
    blue: xr.DataArray,
    *,
    lat_key: str,
    lon_key: str,
    projection: ccrs.Projection,
    title: str,
    output_path: Path,
    plot_options: PlotOptions,
) -> None:
    """Plot a composite RGB map with geographic projection and coastlines.

    Stacks three channels into an RGB image and renders it on a projected
    map axis using imshow with proper geographic extent.

    Args:
        red:          2D DataArray for the red channel (values 0-1).
        green:        2D DataArray for the green channel (values 0-1).
        blue:         2D DataArray for the blue channel (values 0-1).
        lat_key:      Name of the latitude coordinate.
        lon_key:      Name of the longitude coordinate.
        projection:   Cartopy CRS projection for display.
        title:        Plot title string.
        output_path:  Full path for the saved PNG file.
        plot_options: PlotOptions with dpi, figsize, colors, feature toggles.
    """
    facecolor = "black" if plot_options.dark_mode else "white"
    title_color = "white" if plot_options.dark_mode else "black"
    frame_color = "white" if plot_options.dark_mode else "black"

    image = np.dstack((red.values, green.values, blue.values))
    image = np.nan_to_num(image, nan=0.0)
    image = np.clip(image, 0.0, 1.0)

    lons = red.coords[lon_key].values
    lats = red.coords[lat_key].values
    extent = [
        float(lons.min()),
        float(lons.max()),
        float(lats.min()),
        float(lats.max()),
    ]

    fig = plt.figure(figsize=plot_options.figsize, dpi=plot_options.dpi)
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    ax.set_global()
    ax.spines["geo"].set_edgecolor(frame_color)

    _add_map_features(ax, plot_options)

    ax.imshow(
        image,
        origin="lower",
        extent=extent,
        transform=ccrs.PlateCarree(),
        zorder=2,
    )

    ax.set_title(title, color=title_color, fontsize=16)
    fig.savefig(
        output_path,
        dpi=plot_options.dpi,
        bbox_inches="tight",
        facecolor=facecolor,
    )
    plt.close(fig)


# ---------------------------------------------------------------------------
# Data loading functions
# ---------------------------------------------------------------------------


def load_ssha_files(options: Options) -> xr.Dataset:
    """Load SSHA NetCDF files.

    Args:
        options: Options object. Contains:
                     - ssha_folder: Path to SSHA data directory.
                     - tskip:       Load every tskip-th file.

    Returns:
        Combined dataset.
    """
    sshafiles = sorted(options.ssha_folder.glob("*.nc"))
    tskip_files = sshafiles[:: options.tskip]
    logging.info(f"Loading {len(tskip_files)} SSHA files...")

    def _preprocess(ds: xr.Dataset) -> xr.Dataset:
        ds = ds[[options.y_key]]
        if options.xskip > 1:
            ds = ds.isel(
                {
                    options.lat_key: slice(None, None, options.xskip),
                    options.lon_key: slice(None, None, options.xskip),
                }
            )
        ds[options.y_key] = ds[options.y_key].astype(np.float32, copy=False)
        return ds

    input_data = xr.open_mfdataset(
        tskip_files,
        preprocess=_preprocess,
        data_vars="minimal",
        coords="minimal",
        compat="override",
        combine="by_coords",
        parallel=options.mf_parallel,
        chunks={
            options.time_key: options.fake_chunk,
            options.lat_key: options.lat_chunk,
            options.lon_key: options.lon_chunk,
        },
    )
    return input_data


def load_argo_file(options: Options) -> xr.Dataset:
    """Load a single Argo NetCDF file.

    Args:
        options: Options object (argo_folder).

    Returns:
        Argo dataset.
    """
    thefiles = sorted(options.argo_folder.glob("*.nc"))
    thefile = thefiles[0]
    ds = xr.open_dataset(thefile, chunks={})[["ohc_2d", "ohc_2d_700m"]]
    if options.xskip > 1:
        ds = ds.isel(
            {
                options.lat_key: slice(None, None, options.xskip),
                options.lon_key: slice(None, None, options.xskip),
            }
        )
    for v in ("ohc_2d", "ohc_2d_700m"):
        ds[v] = ds[v].astype(np.float32, copy=False)
    return ds


def load_grace_file(options: Options) -> xr.Dataset:
    """Load a single GRACE NetCDF file.

    Args:
        options: Options object (grace_folder).

    Returns:
        GRACE dataset.
    """
    thefiles = sorted(options.grace_folder.glob("*.nc"))
    thefile = thefiles[0]
    ds = xr.open_dataset(thefile, chunks={})[[options.y_key]]
    if options.xskip > 1:
        ds = ds.isel(
            {
                options.lat_key: slice(None, None, options.xskip),
                options.lon_key: slice(None, None, options.xskip),
            }
        )
    ds[options.y_key] = ds[options.y_key].astype(np.float32, copy=False)
    return ds


def load_mur_sst_files(options: Options) -> xr.Dataset:
    """Load MUR SST NetCDF files.

    Args:
        options: Options object (mur_sst_folder, tskip).

    Returns:
        Combined dataset.
    """
    thefiles = sorted(options.mur_sst_folder.glob("*.nc"))
    tskip_files = thefiles[:: options.tskip]
    logging.info(f"Loading {len(tskip_files)} MUR SST files...")

    def _preprocess(ds: xr.Dataset) -> xr.Dataset:
        ds = ds[[options.y_key]]
        if options.xskip > 1:
            ds = ds.isel(
                {
                    options.lat_key: slice(None, None, options.xskip),
                    options.lon_key: slice(None, None, options.xskip),
                }
            )
        ds[options.y_key] = ds[options.y_key].astype(np.float32, copy=False)
        return ds

    input_data = xr.open_mfdataset(
        tskip_files,
        preprocess=_preprocess,
        data_vars="minimal",
        coords="minimal",
        compat="override",
        combine="by_coords",
        parallel=options.mf_parallel,
        chunks={
            options.time_key: options.fake_chunk,
            options.lat_key: options.lat_chunk,
            options.lon_key: options.lon_chunk,
        },
    )
    return input_data


def load_aqua_modis_files(options: Options) -> xr.Dataset:
    """Load AQUA MODIS NetCDF files.

    Args:
        options: Options object (aqua_modis_folder, tskip).

    Returns:
        Combined dataset.
    """
    thefiles = sorted(options.aqua_modis_folder.glob("*.nc"))
    tskip_files = thefiles[:: options.tskip]
    logging.info(f"Loading {len(tskip_files)} AQUA_MODIS files...")

    def _preprocess(ds: xr.Dataset) -> xr.Dataset:
        ds = ds[[options.y_key]]
        if options.xskip > 1:
            ds = ds.isel(
                {
                    options.lat_key: slice(None, None, options.xskip),
                    options.lon_key: slice(None, None, options.xskip),
                }
            )
        ds[options.y_key] = ds[options.y_key].astype(np.float32, copy=False)
        # L3m files carry the date in `time_coverage_start` attr, not a coord.
        time = np.datetime64(ds.attrs["time_coverage_start"][:19])
        ds = ds.assign_coords({options.time_key: time}).expand_dims(options.time_key)
        return ds

    input_data = xr.open_mfdataset(
        tskip_files,
        preprocess=_preprocess,
        combine="nested",
        concat_dim=options.time_key,
        data_vars="minimal",
        coords="minimal",
        compat="override",
        parallel=options.mf_parallel,
        chunks={
            options.time_key: options.fake_chunk,
            options.lat_key: options.lat_chunk,
            options.lon_key: options.lon_chunk,
        },
    )

    return input_data


# Kerchunk per-dataset configuration.  Keys are the dispatch names used by
# _INPUT_DISPATCH; values are MultiZarrToZarr kwargs describing how the
# per-file references should be combined.  Datasets absent from this mapping
# don't have a kerchunk path (single-file datasets, or formats kerchunk
# doesn't yet handle here).
_KERCHUNK_CONFIG: Final[dict[str, dict[str, Any]]] = {
    "AQUA_MODIS": {
        "concat_dim": "time",
        "coo_map": {"time": "attr:time_coverage_start"},
        "coo_dtypes": {"time": "datetime64[ns]"},
        "identical_dims": ["lat", "lon"],
    },
    "SSHA": {"concat_dim": "time"},
    "MUR_SST": {"concat_dim": "time"},
}


def build_kerchunk_reference(
    files: list[Path],
    output_path: Path,
    *,
    concat_dim: str,
    coo_map: dict[str, Any] | None = None,
    coo_dtypes: dict[str, str] | None = None,
    identical_dims: list[str] | None = None,
    inline_threshold: int = 300,
) -> Path:
    """Build a kerchunk reference JSON describing a multi-file NetCDF dataset.

    The reference is a tiny JSON sidecar that lets fsspec's reference
    filesystem expose the on-disk HDF5 chunks as a single virtual Zarr store.
    Once written, ``xr.open_dataset(..., engine="zarr")`` over the reference
    can do chunk-aware random access — which is what allows strided ``isel``
    to reduce real I/O when ``xskip >= disk-chunk-size``.

    Args:
        files: Input NetCDF files in the order they should be concatenated.
        output_path: Destination JSON path.
        concat_dim: Dim along which to concatenate (typically ``"time"``).
        coo_map: kerchunk coo_map for synthesising coords from attrs/data.
        coo_dtypes: kerchunk coo_dtypes for the synthesised coords.
        identical_dims: Dims that must match exactly across files (e.g. lat/lon).
        inline_threshold: Embed array data smaller than this many bytes.

    Returns:
        ``output_path``.
    """
    import json

    # Imported lazily so the regular loader path doesn't pay an import cost.
    from kerchunk.combine import MultiZarrToZarr
    from kerchunk.hdf import SingleHdf5ToZarr

    single_refs = []
    for f in files:
        ref = SingleHdf5ToZarr(str(f), inline_threshold=inline_threshold).translate()
        single_refs.append(ref)

    mzz_kwargs: dict[str, Any] = {"concat_dims": [concat_dim]}
    if coo_map is not None:
        mzz_kwargs["coo_map"] = coo_map
    if coo_dtypes is not None:
        mzz_kwargs["coo_dtypes"] = coo_dtypes
    if identical_dims is not None:
        mzz_kwargs["identical_dims"] = identical_dims

    combined = MultiZarrToZarr(single_refs, **mzz_kwargs).translate()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(combined))
    return output_path


def _kerchunk_reference_path(options: Options) -> Path:
    """Resolve the kerchunk reference file path for the current dataset."""
    if options.kerchunk_path is not None:
        return options.kerchunk_path
    _y_key, folder_attr, _glob = _INPUT_DISPATCH[options.input_choice]
    folder: Path = getattr(options, folder_attr)
    return folder / ".kerchunk.json"


def _kerchunk_reference_is_stale(ref_path: Path, files: list[Path]) -> bool:
    """True iff the reference is missing or older than any input file."""
    if not ref_path.exists():
        return True
    ref_mtime = ref_path.stat().st_mtime_ns
    return any(f.stat().st_mtime_ns > ref_mtime for f in files)


def load_via_kerchunk(options: Options) -> xr.Dataset:
    """Load multi-file NetCDF data via a kerchunk virtual-zarr reference.

    Builds the reference on first call (or when source files are newer) and
    opens the combined dataset through fsspec's reference filesystem.
    The on-disk NetCDF files are not modified; only a small JSON sidecar is
    written next to the data.

    Args:
        options: Options object. Reads ``input_choice``, the matching
                 ``*_folder``, ``kerchunk_path``, ``time_key``, ``lat_key``,
                 ``lon_key``, ``y_key``, ``xskip``, ``lat_chunk``,
                 ``lon_chunk``, ``fake_chunk``, ``tskip``.

    Returns:
        Combined dataset shaped like the matching non-kerchunk loader would
        have produced (after preprocess select / xskip / cast).

    Raises:
        ValueError: If ``input_choice`` has no kerchunk configuration.
        FileNotFoundError: If no input files match the dispatch glob.
    """
    import fsspec

    if options.input_choice not in _KERCHUNK_CONFIG:
        raise ValueError(
            f"input_choice {options.input_choice!r} has no kerchunk config; "
            f"use the regular loader."
        )

    _y_key, folder_attr, glob_pattern = _INPUT_DISPATCH[options.input_choice]
    folder: Path = getattr(options, folder_attr)
    all_files = sorted(folder.glob(glob_pattern))
    if not all_files:
        raise FileNotFoundError(f"No files matching {glob_pattern!r} in {folder}")
    tskip_files = all_files[:: options.tskip]

    ref_path = _kerchunk_reference_path(options)
    if _kerchunk_reference_is_stale(ref_path, tskip_files):
        cfg = _KERCHUNK_CONFIG[options.input_choice]
        logging.info(
            f"Building kerchunk reference for {len(tskip_files)} "
            f"{options.input_choice} files → {os.fspath(ref_path)}"
        )
        build_kerchunk_reference(tskip_files, ref_path, **cfg)
    else:
        logging.info(f"Reusing kerchunk reference at {os.fspath(ref_path)}")

    fs = fsspec.filesystem("reference", fo=str(ref_path))
    mapper = fs.get_mapper("")
    ds = xr.open_dataset(
        mapper,
        engine="zarr",
        consolidated=False,
        chunks={
            options.time_key: options.fake_chunk,
            options.lat_key: options.lat_chunk,
            options.lon_key: options.lon_chunk,
        },
    )

    ds = ds[[options.y_key]]
    if options.xskip > 1:
        ds = ds.isel(
            {
                options.lat_key: slice(None, None, options.xskip),
                options.lon_key: slice(None, None, options.xskip),
            }
        )
    ds[options.y_key] = ds[options.y_key].astype(np.float32, copy=False)
    return ds


def load_simple_grid_files(options: Options) -> xr.Dataset:
    """Load simple_grids weekly NetCDF files and concatenate along time.

    Each file contains a single weekly snapshot of ssha(latitude, longitude)
    with a scalar 'time' data variable.  Files are found recursively under
    options.simple_folder (which contains year subfolders like 2025/, 2026/).

    Args:
        options: Options object. Contains:
                     - simple_folder: Path to simple_grids data directory.
                     - tskip:         Load every tskip-th file (temporal decimation).

    Returns:
        Combined dataset with dimensions (time, latitude, longitude).

    Raises:
        FileNotFoundError: If no NetCDF files are found in simple_folder.
    """
    all_files = sorted(options.simple_folder.rglob("*.nc"))
    if not all_files:
        raise FileNotFoundError(
            f"No NetCDF files found in {os.fspath(options.simple_folder)}"
        )
    tskip_files = all_files[:: options.tskip]
    logging.info(
        f"Loading {len(tskip_files)} simple_grids files (of {len(all_files)} total)..."
    )

    def _preprocess(ds: xr.Dataset) -> xr.Dataset:
        """Promote scalar time variable to a dimension coordinate."""
        keep = {options.y_key, options.time_key}
        drop = set(ds.data_vars) - keep
        ds = ds.drop_vars(drop, errors="ignore")
        if options.xskip > 1:
            ds = ds.isel(
                {
                    options.lat_key: slice(None, None, options.xskip),
                    options.lon_key: slice(None, None, options.xskip),
                }
            )
        ds = ds.set_coords(options.time_key).expand_dims(options.time_key)
        ds[options.y_key] = ds[options.y_key].astype(np.float32, copy=False)
        return ds

    input_data = xr.open_mfdataset(
        tskip_files,
        preprocess=_preprocess,
        combine="nested",
        concat_dim=options.time_key,
        data_vars="minimal",
        coords="minimal",
        compat="override",
        parallel=options.mf_parallel,
        chunks={
            options.time_key: options.fake_chunk,
            options.lat_key: options.lat_chunk,
            options.lon_key: options.lon_chunk,
        },
    )
    return input_data


# ---------------------------------------------------------------------------
# GMT high-level functions
# ---------------------------------------------------------------------------


def write_gmt_cpt(options: Options, new_fp: BinaryIO, i: int) -> None:
    """Write GMT commands for creating the CPT (color palette) file.

    Ported from C++ write_gmt_cpt() (abridged-functions.cpp lines 3-98).
    Only called for non-RGB maps (single-parameter: trend, amplitude, phase).
    For RGB spectral-color maps, this function is not invoked.

    Args:
        options: Options object. Contains:
                     - results:      Results with options.output_choice, rgb, etc.
                     - plot_options: PlotOptions with symmetric_limit, scale_digits.
        new_fp:  Binary file handle for the GMT script being written.
        i:       Index of the parameter being scaled.
    """
    results = options.results
    plot_options = options.plot_options
    numlevels = 10
    digits = plot_options.scale_digits

    palette_guide = (
        "#0-360phase=cyclic,amps=wysiwyg,trend=polar,topo=relief,pts=rainbow"
    )
    overflow_guide = "#-E=both,-Ef=top,-Eb=bottom"

    if results.options.output_choice == 1:
        # Unstructured output
        if results.options.output_type[i] == 104:
            # Phase map: cyclic palette, fixed [-180, 180] range
            new_fp.write(f'palette=" -Ccyclic " {palette_guide}\n'.encode())
            new_fp.write(f'upper_limit="{180.0:.6f}"\n'.encode())
            new_fp.write(f'lower_limit="{-180.0:.6f}"\n'.encode())
        else:
            new_fp.write(f'palette=" -Crainbow " {palette_guide}\n'.encode())
            new_fp.write(f'upper_limit="{results.max:.6f}"\n'.encode())
            new_fp.write(f'lower_limit="{results.min:.6f}"\n'.encode())
        new_fp.write(f'overflow="" {overflow_guide}\n'.encode())
        new_fp.write(f'numlevels="{numlevels}"\n'.encode())
        new_fp.write(b'# The "| bc -l" helps bash deal with floating point numbers.\n')
        new_fp.write(
            b'interval="$(echo "($upper_limit - $lower_limit) / $numlevels" | bc -l)"\n'
        )
        new_fp.write(
            b"$gmt_prefix makecpt $palette -T$lower_limit/$upper_limit/$interval -Z > map.cpt\n"
        )

    elif results.options.output_choice == 5:
        # Latlon grid output — check RGB first since output_type may be empty
        if len(results.rgb) == 3:
            # RGB map: just set coast color
            new_fp.write(b'  coast_color="gray82"\n')

        elif (
            i < len(results.options.output_type)
            and results.options.output_type[i] == 104
        ):
            # Phase map: cyclic palette
            logging.info(
                "NOTE: Phases should only use the cyclic colorscale if the min/max span [-180,180] or [0,360]!"
            )
            new_fp.write(f'palette=" -Ccyclic -I " {palette_guide}\n'.encode())
            new_fp.write(
                b'limits="" #Default limits include all values (not just in subset), aren\'t always symmetric.\n'
            )
            new_fp.write(f'overflow="" {overflow_guide}\n'.encode())
            numlevels = 13  # 360 degrees / 12 (+1 for value 0.0) divides nicely
            new_fp.write(f'numlevels="{numlevels}"\n'.encode())
            new_fp.write(
                b"$gmt_prefix grd2cpt $data_name $palette $limits -E$numlevels -Z > map.cpt\n"
            )

        else:
            # Non-RGB single-parameter map: WYSIWYG palette with data-driven limits
            new_fp.write(
                b'subset_name=$data_name"_subset" #Weird to append, but works with ../pl..\n'
            )
            new_fp.write(
                b"$gmt_prefix grdcut $data_name -G$subset_name $actual_range #Only use subset values for scale.\n"
            )
            new_fp.write(
                b"$gmt_prefix grdinfo -L0 $subset_name > subset_grdinfo #Extract data range from subset.\n"
            )
            new_fp.write(
                f"data_min_e=$(awk '/z_min: /{{printf \"%.{digits}e\\n\", $3}}' subset_grdinfo)\n".encode()
            )
            new_fp.write(
                f"data_max_e=$(awk '/z_max: /{{printf \"%.{digits}e\\n\", $5}}' subset_grdinfo)\n".encode()
            )
            new_fp.write(
                f"data_min_f=$(awk '/z_min: /{{printf \"%.{digits}f\\n\", $3}}' subset_grdinfo)\n".encode()
            )
            new_fp.write(
                f"data_max_f=$(awk '/z_max: /{{printf \"%.{digits}f\\n\", $5}}' subset_grdinfo)\n".encode()
            )
            new_fp.write(
                b'notation="f" #By default, numbers appear in floating point format.\n'
            )
            new_fp.write(
                b"data_min_print=$data_min_f #By default, numbers appear in floating point format.\n"
            )
            new_fp.write(
                b"data_max_print=$data_max_f #By default, numbers appear in floating point format.\n"
            )
            new_fp.write(f'palette=" -Cwysiwyg " {palette_guide}\n'.encode())
            new_fp.write(
                f'symmetric_limit="{plot_options.symmetric_limit:.6f}"\n'.encode()
            )
            new_fp.write(
                b'if [[ $(bc <<< "$symmetric_limit <= 0.0") == 1 || $(bc <<< "$data_min_f >= 0.0") == 1 || $(bc <<< "$data_max_f <= 0.0 ") == 1 ]]\n'
            )
            new_fp.write(b"then\n")
            new_fp.write(b"  upper_limit=$data_max_f\n")
            new_fp.write(b"  lower_limit=$data_min_f\n")
            new_fp.write(b"else\n")
            new_fp.write(b"  upper_limit=$symmetric_limit\n")
            new_fp.write(b"  lower_limit=-$upper_limit\n")
            new_fp.write(b"fi\n")
            new_fp.write(b'limits=" -L$lower_limit/$upper_limit "\n')
            new_fp.write(b". ./overflow.sh\n")
            # Coastline color: off-white if all non-negative, dark if mixed signs
            new_fp.write(b'if [[ $(bc <<< "$data_min_f >= 0.0") == 1 ]]\n')
            new_fp.write(b"then\n")
            new_fp.write(b'  coast_color="gray82"\n')
            new_fp.write(b"else\n")
            new_fp.write(b'  coast_color="gray10"\n')
            new_fp.write(b"fi\n")
            new_fp.write(f'numlevels="{numlevels}"\n'.encode())
            new_fp.write(
                b"$gmt_prefix grd2cpt $data_name $palette $limits -E$numlevels -Z > map.cpt\n"
            )

    else:
        logging.error(
            f"!!!!WARNING!!!!!! results.options.output_choice {results.options.output_choice} isn't recognized."
        )

    plot_options.blurb_disabled = 0  # Better to always see the blurb


def is_polar(options: Options) -> int:
    """Determine if the data is global, north-polar, or south-polar.

    Args:
        options: Options object. Contains:
                     - results: Results (mutated with min/maxlat/lon).
                     - grid:    GridData.

    Returns:
        0 for global, 1 for north pole, 2 for south pole.
    """
    results = options.results
    grid = options.grid
    polar = 0

    if results.options.output_choice in (1, 4):
        results.minlat = min(grid.lat)
        results.maxlat = max(grid.lat)
        results.minlon = min(grid.lon)
        results.maxlon = max(grid.lon)
    elif results.options.output_choice == 5:
        results.minlat = min(results.latlon.lat)
        results.maxlat = max(results.latlon.lat)
        results.minlon = min(results.latlon.lon)
        results.maxlon = max(results.latlon.lon)
    else:
        logging.error(
            f"!!!!WARNING!!!!!! results.options.output_choice {results.options.output_choice} isn't recognized."
        )

    if results.minlat < 0 and results.maxlat > 0:
        polar = 0
    elif results.minlat > 0 and results.maxlat > 0:
        polar = 1
    elif results.minlat < 0 and results.maxlat < 0:
        polar = 2

    return polar


def write_gmt_colorscale(options: Options, new_fp: BinaryIO, kml_output: int) -> None:
    """Write the GMT color scale command to the script file.

    Args:
        options:    Options object (results).
        new_fp:     File handle to the GMT script.
        kml_output: 0 for GMT scale next to plot, 1 for KMZ by itself.
    """
    results = options.results

    if len(results.rgb) != 3:
        new_fp.write(b'cpt_name="-Cmap.cpt "\n')
    else:
        new_fp.write(b'cpt_name="-Crgb00001.cpt "\n')

    if len(results.rgb) != 3 or results.rgb_choice < 2:
        if kml_output:
            new_fp.write(
                b"gmt psscale $cpt_name -L $scale_format $overflow $scale_pos -A $start > scale_$plot_base.ps\n"
            )
            new_fp.write(
                b"# Print units manually, otherwise they're too close to numbers on scale.\n"
            )
            new_fp.write(
                b"echo $units_format $scale_units | gmt pstext -N $units_pos $misc_range $end >> scale_$plot_base.ps\n"
            )
        else:
            new_fp.write(
                b"gmt psscale $cpt_name -L $scale_format $overflow $scale_pos -A $middle >> $plot_base.ps\n"
            )
            new_fp.write(
                b"# Print units manually, otherwise they're too close to numbers on scale.\n"
            )
            new_fp.write(
                b"echo $units_format $scale_units | gmt pstext -N $units_pos $misc_range $middle >> $plot_base.ps\n"
            )
    else:
        new_fp.write(f"numwidths={results.max_widths}\n".encode())
        new_fp.write(b"gmt set TICK_LENGTH 0.3c\n")
        new_fp.write(b'scale_width=$(bc <<< "scale=5; $scale_width / $numwidths")\n')
        new_fp.write(b"for (( j=1; j <= $numwidths; j++ ))\n")
        new_fp.write(b"do\n")
        new_fp.write(b"  j_string=$(printf '%%05d' $j)\n")
        new_fp.write(b'  cpt_name="-Crgb$j_string.cpt"\n')
        new_fp.write(
            b"  gmt psscale $cpt_name -L $scale_format $overflow $scale_pos -S -A $middle >> $plot_base.ps\n"
        )
        new_fp.write(
            b"  # Print units manually, otherwise they're too close to numbers on scale.\n"
        )
        new_fp.write(
            b"  echo $units_format $scale_units | gmt pstext -N $units_pos $misc_range $middle >> $plot_base.ps\n"
        )
        new_fp.write(
            b"  # Move next scale to the right and get rid of the tick marks.\n"
        )
        new_fp.write(b'  scale_x=$(bc <<< "scale=5; $scale_x+$scale_width")\n')
        new_fp.write(
            b'  scale_pos=" -D${scale_x}c/${scale_y}c/${scale_length}c/${scale_width}c "\n'
        )
        new_fp.write(b"  gmt set TICK_LENGTH 0.0\n")
        new_fp.write(b"done\n")


def write_rgb_colorscale(options: Options, cie: xr.Dataset) -> None:
    """Generate a GMT CPT file mapping period values to spectral RGB colors.

    For each evenly-spaced period in [min_period, max_period], computes the
    spectral color that a pure sinusoid at that period would produce using
    the same CIE color-matching pipeline as the main hot path.  Writes the
    result as a GMT-format CPT file (rgb00001.cpt).

    Args:
        options: Options object. Contains:
                     - min_period:    Minimum period in days.
                     - max_period:    Maximum period in days.
                     - thepower:      Y-brightening exponent.
                     - use_new_funcs: True for colour-science, False for Riemann sum.
                     - output_folder: Output directory Path.
        cie:     CIE 1931 color matching functions Dataset with 'x', 'y', 'z'
                 variables and 'wavelength' coordinate.
    """
    n_steps = 256
    cie_wavelengths = cie["wavelength"].values
    wl_min = cie_wavelengths[0]
    wl_range = cie_wavelengths[-1] - wl_min

    # Sample periods and their corresponding CIE wavelengths
    periods = np.linspace(options.min_period, options.max_period, n_steps)
    target_wavelengths = (
        wl_min
        + (periods - options.min_period)
        / (options.max_period - options.min_period)
        * wl_range
    )

    # Build power matrix: each row is a narrow Gaussian centered at the target wavelength
    sigma = 1.0  # nm — narrow enough to approximate a delta function
    power_matrix = np.exp(
        -np.power(
            cie_wavelengths[np.newaxis, :] - target_wavelengths[:, np.newaxis], 2.0
        )
        / (2.0 * sigma**2)
    )

    # Compute XYZ for each sample period
    xyz_cpt = np.zeros((n_steps, 3))

    if options.use_new_funcs:
        cmfs = MSDS_CMFS_STANDARD_OBSERVER["CIE 1931 2 Degree Standard Observer"]
        for i in range(n_steps):
            spd = SpectralDistribution(power_matrix[i], cie_wavelengths)
            xyz_cpt[i] = sd_to_XYZ(spd, cmfs)
    else:
        cie_x_vals = cie["x"].values
        cie_y_vals = cie["y"].values
        cie_z_vals = cie["z"].values
        wavelength_step = np.diff(cie_wavelengths).mean()
        xyz_cpt[:, 0] = (power_matrix * cie_x_vals[np.newaxis, :]).sum(
            axis=1
        ) * wavelength_step
        xyz_cpt[:, 1] = (power_matrix * cie_y_vals[np.newaxis, :]).sum(
            axis=1
        ) * wavelength_step
        xyz_cpt[:, 2] = (power_matrix * cie_z_vals[np.newaxis, :]).sum(
            axis=1
        ) * wavelength_step

    # --- Post-processing (same as hot path lines 338-373) ---

    # Normalize by max Y
    max_y = xyz_cpt[:, 1].max()
    xyz_cpt /= max_y

    # Raise Y to power while keeping chromaticity constant
    factor = np.power(xyz_cpt[:, 1], 1.0 - options.thepower)
    xyz_cpt /= factor[:, np.newaxis]

    # XYZ to sRGB
    if options.use_new_funcs:
        rgb_cpt = XYZ_to_sRGB(xyz_cpt)
    else:
        rgb_cpt = (A_OLD_XYZ_TO_SRGB @ xyz_cpt.T).T

    # Fix gamut (vectorized)
    rgb_arr = rgb_cpt.T.copy()  # shape (3, n_steps)
    min_vals = rgb_arr.min(axis=0)
    needs_fix = min_vals < 0

    if np.any(needs_fix):
        offset = -min_vals[needs_fix] + 1.0 / 255.0
        y_vals = xyz_cpt[needs_fix, 1]
        fix_factor = y_vals / (y_vals + offset)
        rgb_arr[:, needs_fix] = (rgb_arr[:, needs_fix] + offset) * fix_factor

    # Normalize to [0, 1]
    rgb_arr /= np.max(rgb_arr)

    # Apply BT.709 gamma correction (Hughes & Williams 2010, eqs. A7-A8)
    if not options.use_new_funcs:
        gamma_inv = 0.45
        crit = 0.018
        h = 4.506813168
        g_coeff = -0.09914989
        f_coeff = 1.09914989
        for ch in range(3):
            vals = rgb_arr[ch]
            low = vals <= crit
            high = vals > crit
            vals[low] = vals[low] * h
            vals[high] = f_coeff * np.power(vals[high], gamma_inv) + g_coeff

    # Scale to [0, 255]
    rgb_arr *= 255.0
    rgb_arr = np.clip(rgb_arr, 0.0, 255.0).astype(int)

    # Write CPT file
    cpt_path = options.output_folder / "rgb00001.cpt"
    logging.info(f"Writing RGB colorscale to {os.fspath(cpt_path)}")
    with cpt_path.open("w") as fp:
        fp.write("# COLOR_MODEL = RGB\n")
        for i in range(1, n_steps):
            fp.write(
                f"{periods[i - 1]:12.6e} {rgb_arr[0, i - 1]:3d} {rgb_arr[1, i - 1]:3d} {rgb_arr[2, i - 1]:3d} "
                f"{periods[i]:12.6e} {rgb_arr[0, i]:3d} {rgb_arr[1, i]:3d} {rgb_arr[2, i]:3d}\n"
            )

    logging.info("Finished writing RGB colorscale.")


def write_gmt_map_data(options: Options, new_fp: BinaryIO, title: str, i: int) -> None:
    """Write GMT data plotting commands to the script file.

    Args:
        options: Options object (results, grid, plot_options).
        new_fp:  File handle to the GMT script.
        title:   Map title string. If blank, outputs KMZ format.
        i:       Index of the parameter being mapped.
    """
    results = options.results
    plot_options = options.plot_options

    polar = is_polar(options)
    coastlines = plot_options.coastlines  # Default; overridden below for polar == 1

    if i == 0:
        new_fp.write(b"#Set resolution, coast_file, coast_thickness, and coastlines\n")
        new_fp.write(b"#on first map only because they should be universal.\n")
        coast_file = (
            plot_options.outputfolder + "data/ancillary/Rignot/InSAR_GL_Antarctica.txt"
        )
        new_fp.write(f'coast_file="{coast_file}"\n'.encode())

        if results.options.output_choice == 5:
            delta_lat = abs(results.latlon.lat[1] - results.latlon.lat[0])
            if delta_lat < 0.4:
                new_fp.write(b'resolution=" -E50 " #50/2000 is low/high quality.\n')
                new_fp.write(b'coast_res=" -Df+ "\n')
            elif delta_lat < 0.9:
                new_fp.write(b'resolution=" -E50 " #50/2000 is low/high quality.\n')
                new_fp.write(b'coast_res=" -Df+ "\n')
            else:
                new_fp.write(b'resolution=" -E50 " #50/2000 is low/high quality.\n')
                new_fp.write(b'coast_res=" -Di+ "\n')
        elif results.options.output_choice in (1, 4):
            new_fp.write(b'resolution=" -E50 " #50/2000 is low/high quality.\n')
            new_fp.write(b'coast_res=" -Di+ "\n')
        else:
            logging.error(
                f"!!!!WARNING!!!!!! results.options.output_choice {results.options.output_choice} isn't recognized."
            )

        new_fp.write(
            b"coast_res_orig=$coast_res #Don't want USA maps to repeatedly add -N2.\n"
        )
        new_fp.write(b'coast_thk="0.6"\n')
        new_fp.write(b'coast_thk="0.009"\n')
        if polar == 1:
            coastlines = 1  # InSAR is only in Antarctica, so disable for NP plots.
        new_fp.write(f"coastlines={coastlines} #1:coast, 2:coast+InSAR.\n".encode())

    new_fp.write(
        b"#coast_color is gray82 for off-white, or gray10 for dark coastlines.\n"
    )

    # Adjust max/min latitudes for mapping points
    if results.options.output_choice in (1, 4):
        buffer = 5
        if results.maxlat <= 90 - buffer:
            results.maxlat += buffer
        if results.minlat >= -90 + buffer:
            results.minlat -= buffer

    new_fp.write(b'title_format="0 0 30 0 0 MC"\n')
    new_fp.write(b'blurb_format="0 0 15 0 1 ML"\n')
    new_fp.write(b'units_format="0 0 13 0 0 MC"\n')
    new_fp.write(f'scale_units="{results.units[i]}"\n'.encode())

    new_fp.write(b'misc_range=" -R0/1/0/1 -JX1c "\n')
    new_fp.write(
        b"#grdcut requires actual limits, but if grdimage uses them: GMT Fatal Error: grdimage could not allocate memory [21.69 Gb, n_items = 5823567396]\n"
    )
    new_fp.write(b"minlon=%.3f\n" % 0.0)
    new_fp.write(b"maxlon=%.3f\n" % 360.0)
    new_fp.write(b"minlat=%.3f\n" % -90.0)
    new_fp.write(b"maxlat=%.3f\n" % 90.0)

    if title:
        if i == 0:
            new_fp.write(b"#Global projections:\n")
            new_fp.write(b"#    1 - Robinson\n")
            new_fp.write(b"#    2 - Winkel Tripel\n")
            new_fp.write(b"#    3 - Mollweide\n")
            new_fp.write(b"#    4 - Miller\n")
            new_fp.write(b"#Polar projections:\n")
            new_fp.write(b"#  101 - N. Azimuthal Equidistant\n")
            new_fp.write(b"#  102 - S. Azimuthal Equidistant\n")
            new_fp.write(b"#Specific regions:\n")
            new_fp.write(b"# 1001 - North America\n")
            new_fp.write(b"# 1002 - South America\n")
            new_fp.write(b"# 1003 - Africa\n")
            new_fp.write(b"# 1004 - Greenland\n")
            new_fp.write(b"# 1005 - South Asia\n")
            new_fp.write(b"# 1006 - Australia\n")
            new_fp.write(b"# 1007 - Europe\n")
            new_fp.write(b"# 1101 - Contiguous United States\n")
            new_fp.write(b"# 1102 - California\n")

        if polar == 0:
            new_fp.write(
                f"{'#' if i > 0 else ''}projection_choice={plot_options.projection}\n".encode()
            )
        else:
            if polar == 1:
                new_fp.write(f"{'#' if i > 0 else ''}projection_choice=101\n".encode())
            else:
                new_fp.write(f"{'#' if i > 0 else ''}projection_choice=102\n".encode())

        new_fp.write(
            b"standard_circle=0 #1=all specific regions use standard circular projection.\n"
        )
        new_fp.write(
            b"standard_rect=0 #1=all specific regions use standard rectangular projection.\n"
        )
        new_fp.write(b". ./projections.sh\n")
        new_fp.write(b"if [ $projection_choice == 101 ]\n")
        new_fp.write(b"then\n")

        minlat = results.minlat if polar == 1 else 0.0
        new_fp.write(f"  minlat={minlat:.3f}\n".encode())
        new_fp.write(b'  actual_range=" -R0.0/360.0/$minlat/90.0 "\n')
        polar_radius = 90 - minlat
        new_fp.write(f"  polar_radius={polar_radius}\n".encode())
        new_fp.write(
            b'  projection=" -JE0/90.0/${polar_radius}/${map_width}c " #N. Azimuthal Equidistant\n'
        )
        new_fp.write(b"elif [ $projection_choice == 102 ]\n")
        new_fp.write(b"then\n")

        maxlat = results.maxlat if polar == 2 else 0.0
        new_fp.write(f"  maxlat={maxlat:.3f}\n".encode())
        new_fp.write(b'  actual_range=" -R0.0/360.0/-90.0/$maxlat "\n')
        polar_radius = 90 + maxlat
        new_fp.write(f"  polar_radius={polar_radius}\n".encode())
        new_fp.write(
            b'  projection=" -JE0/-90.0/{polar_radius}/${map_width}c " #S. Azimuthal Equidistant\n'
        )
        new_fp.write(b"fi\n")

        new_fp.write(b'range=" -R${minlon}/${maxlon}/${minlat}/${maxlat} "\n')
        new_fp.write(b'map_pos=" -Xa${map_x}c -Ya${map_y}c "\n')

        if len(results.rgb) == 3 and results.rgb_choice >= 2:
            new_fp.write(b"scale_width=1.2 #Override for RGB maps.\n")

        new_fp.write(
            b'scale_pos=" -D${scale_x}c/${scale_y}c/${scale_length}c/${scale_width}c "\n'
        )
        new_fp.write(b'units_x=$(bc <<< "scale=5; $scale_x+$scale_width/2")\n')
        new_fp.write(b'units_y=$(bc <<< "scale=5; $scale_y+$scale_length/2")\n')
        new_fp.write(b'units_pos=" -Xa${units_x}c -Ya${units_y}c "\n')
        new_fp.write(b'blurb_pos=" -Xa${blurb_x}c -Ya${blurbs_y}c "\n')
        new_fp.write(b'blurb2_pos=" -Xa${blurb2_x}c -Ya${blurbs_y}c "\n')
    else:
        logging.error("!!!WARNING!!! NO TITLE!")

    # Create CPT file for coloring map (non-RGB maps only).
    write_gmt_cpt(options, new_fp, i)

    # Plot data, with title on top.
    new_fp.write(f'title="{title}"\n'.encode())
    if len(results.rgb) == 3 and len(results.latlon.outputs) == 3:
        new_fp.write(
            b"gmt grdinfo red_gmt.nc?z || echo 'ERROR: GMT cannot read red_gmt.nc' >&2\n"
        )
        new_fp.write(
            b"gmt grdimage red_gmt.nc?z green_gmt.nc?z blue_gmt.nc?z $boundary $resolution $range $projection $map_pos $start > $plot_base.ps\n"
        )
        new_fp.write(
            b'if [ $? -ne 0 ]; then echo "ERROR: gmt grdimage failed" >&2; fi\n'
        )
    else:
        new_fp.write(
            b"gmt grdimage $data_name $boundary $resolution $range $projection $map_pos -Cmap.cpt $start > $plot_base.ps\n"
        )

    write_gmt_coastlines(new_fp)

    if (
        len(results.marker_lats) == len(results.latlon.outputs)
        and len(results.marker_lons) == len(results.latlon.outputs)
        and results.latlon.outputs
    ):
        if results.marker_lats[i] and results.marker_lons[i]:
            new_fp.write(
                b"gmt psxy -N $data_name -bcmarker_lons/marker_lats -S+0.5c -W5/244/164/96 -G244/164/96 $range $projection $map_pos $middle >> $plot_base.ps\n"
            )

    new_fp.write(b"#Uncomment to put a marker at echoed coords, given as lon lat:\n")
    new_fp.write(
        b"#echo -85.19 -77.36 | gmt psxy -N -S+0.5c -W5/244/164/96 -G244/164/96 $range $projection $map_pos $middle >> $plot_base.ps\n"
    )

    if plot_options.plot_mascons != 0 and results.latlon.mascon_lats:
        new_fp.write(
            b"gmt psxy $data_name -bcmascon_lons/mascon_lats -Sc0.01c -G139/69/19 $range $projection $map_pos $middle >> $plot_base.ps\n"
        )

    new_fp.write(b"if [ $montage != 0 ]\n")
    new_fp.write(b"then\n")
    new_fp.write(b'  title=${prefixes[$index]}" "$title\n')
    new_fp.write(
        b'  title_format="0 0 30 0 0 ML" #Left-justify so montage titles are uniform.\n'
    )
    new_fp.write(b'  title_x=$(bc <<< "scale=5; $blurb_x-0.1")\n')
    new_fp.write(b"else\n")
    new_fp.write(b'  title_x=$(bc <<< "scale=5; $map_x+$map_width/2")\n')
    new_fp.write(b"fi\n")
    new_fp.write(b'title_pos=" -Xa${title_x}c -Ya${title_y}c "\n')
    new_fp.write(
        b"echo $title_format $title | gmt pstext -N $title_pos $misc_range $middle >> $plot_base.ps\n"
    )

    # Draw color scale with units printed above.
    write_gmt_colorscale(options, new_fp, 0)

    # Print blurb about data range, or masked amplitudes for phase plots.
    if len(results.rgb) == 3 and len(results.latlon.outputs) == 3:
        plot_options.blurb_disabled = 1
    if plot_options.blurb_disabled:
        new_fp.write(b'blurb_contents=""\n')

    new_fp.write(
        b"echo $blurb_format $blurb_contents | gmt pstext -N $blurb_pos $misc_range $middle >> $plot_base.ps\n"
    )

    blurb2_written = 0
    if (
        len(results.latlon.outputs) == len(results.error_bars)
        and results.latlon.outputs
    ):
        if results.error_bars[i] > 0.0:
            blurb2_written = 1
            new_fp.write(
                f'blurb2_contents="Error bar: {results.error_bars[i]:.1f} $scale_units"\n'.encode()
            )

    if not blurb2_written:
        new_fp.write(b'blurb2_contents="" #Error bar: N/A $scale_units\n')

    new_fp.write(
        b"echo $blurb_format $blurb2_contents | gmt pstext -N $blurb2_pos $misc_range $end >> $plot_base.ps\n"
    )


def write_gmt_scripts(options: Options) -> None:
    """Generate GMT shell scripts for creating map plots.

    Args:
        options: Options object (plot_options, grid, results).
    """
    plot_options = options.plot_options
    results = options.results

    new_file = Path(plot_options.outputfolder) / "create_plots.sh"
    try:
        new_fp = new_file.open("wb")
    except OSError:
        logging.error("The create_plots.sh GMT script couldn't be created.")

    new_fp.write(b"#!/bin/bash\n")
    new_fp.write(b"#set -x #Uncomment to echo these commands.\n")

    flip_file = Path(plot_options.outputfolder) / "flip_backgrounds.sh"
    try:
        flip_fp = flip_file.open("wb")
    except OSError:
        logging.error("The flip_backgrounds.sh script couldn't be created.")

    flip_fp.write(b"#!/bin/bash\n")
    flip_fp.write(b"set -x\n")

    trim_file = Path(plot_options.outputfolder) / "trim.sh"
    try:
        trim_fp = trim_file.open("wb")
    except OSError:
        logging.error("The trim.sh script couldn't be created.")

    trim_fp.write(b"#!/bin/bash\n")
    trim_fp.write(b"set -x\n")

    if len(results.rgb) == 3 and len(results.latlon.outputs) == 3:
        just_the_filenames = plot_options.just_the_filenames[:1]
    else:
        just_the_filenames = plot_options.just_the_filenames

    for i in range(len(just_the_filenames)):
        if i == 0:
            write_gmt_defs(new_fp)
            new_fp.write(
                f"color_scheme={plot_options.color_scheme} #1/2=white/black background\n".encode()
            )
            new_fp.write(
                f"montage={plot_options.montage} #1=left-justify titles, add (a),(b), run montage.sh.\n".encode()
            )
            new_fp.write(
                b'prefixes=("(a)" "(b)" "(c)" "(d)" "(e)" "(f)" "(g)" "(h)" "(i)" "(j)" "(k)" "(l)" "(m)" "(n)" "(o)" "(p)" "(q)" "(r)" "(s)" "(t)" "(u)" "(v)" "(w)" "(x)" "(y)" "(z)")\n'
            )
            new_fp.write(
                b"index=-1 #Increments on each map, accesses prefixes above for montage.\n"
            )
            new_fp.write(b'png_options=" -P -Tg " #PDF default: -E720, else 300 dpi.\n')
            new_fp.write(
                b'if [ $montage != 0 ]\nthen\n  png_options=" -A"$png_options\nfi\n'
            )
            new_fp.write(
                b"#Force off-white(dark gray) fore(back)ground color because\n#flip_backgrounds.sh can change the maps' text from\n#black to white, and their backgrounds from white to black.\n"
            )
            new_fp.write(
                b"gmt set COLOR_BACKGROUND=2/2/2 COLOR_FOREGROUND=253/253/253\n"
            )
            new_fp.write(f"digits={plot_options.scale_digits}\n".encode())
            new_fp.write(b"gmt set D_FORMAT=%.${digits}f\n")

        s = f"{plot_options.output_base}_{i + 1:04d}"
        new_fp.write(b"#######################################################\n")
        if len(results.rgb) == 3 and len(results.latlon.outputs) == 3:
            new_fp.write(b"data_name=redgreenblue\n")
        else:
            new_fp.write(f'data_name="{just_the_filenames[i]}"\n'.encode())
        new_fp.write(f'plot_base="{s}"\n'.encode())
        new_fp.write(b"let index=$index+1\n")
        new_fp.write(b"#######################################################\n")

        if i == 0:
            if len(just_the_filenames) == 1:
                flip_fp.write(f'all_bases="{s}"\n'.encode())
                trim_fp.write(f'all_bases="{s}"\n'.encode())
            else:
                flip_fp.write(f'all_bases="{s}\n'.encode())
                trim_fp.write(f'all_bases="{s}\n'.encode())
        elif i < len(just_the_filenames) - 1:
            flip_fp.write(f"{s}\n".encode())
            trim_fp.write(f"{s}\n".encode())
        else:
            flip_fp.write(f'{s}"\n'.encode())
            trim_fp.write(f'{s}"\n'.encode())

        write_gmt_map_data(options, new_fp, results.titles[i], i)

        new_fp.write(
            b"gmt psconvert $png_options $plot_base.ps #Convert PS to PNG format.\n"
        )
        new_fp.write(
            b"#convert -P -Tf $plot_base.ps #Convert PS to PDF, if uncommented.\n"
        )
        new_fp.write(b"rm -f $plot_base.ps\n")

        if len(results.rgb) != 3:
            new_fp.write(b"mv map.cpt Zbackup_cpt_$plot_base.cpt\n")

        new_fp.write(b"previous_data_name=$data_name\n")

    new_fp.write(b"#######################################################\n")
    new_fp.write(b". ./trim.sh\n")

    new_fp.write(b"#######################################################\n")
    new_fp.write(b"if [ $color_scheme == 2 ]\n")
    new_fp.write(b"then\n")
    new_fp.write(b"  . ./flip_backgrounds.sh\n")
    new_fp.write(b"fi\n")

    new_fp.write(b"#######################################################\n")
    new_fp.write(b"if [ $montage != 0 ]\n")
    new_fp.write(b"then\n")
    new_fp.write(b"  . ./montage.sh\n")
    new_fp.write(b"fi\n")
    new_fp.close()

    finish_flip_backgrounds(flip_fp)
    flip_fp.close()

    finish_trim(trim_fp)
    trim_fp.close()

    # Create animate script
    extra_file = Path(plot_options.outputfolder) / "animate.sh"
    try:
        extra_fp = extra_file.open("wb")
    except OSError:
        logging.error("The animate.sh script couldn't be created.")

    extra_fp.write(b"#!/bin/bash\n")
    extra_fp.write(b"set -x\n")
    extra_fp.write(b"delay=100 #delay in hundredths of a second.\n")
    extra_fp.write(b'#size="640x480"\n')
    extra_fp.write(b'#size="800x600"\n')
    extra_fp.write(b'size="1024x768"\n')
    extra_fp.write(f'output_base="{plot_options.output_base}"\n'.encode())
    extra_fp.write(
        b"#Imagemagick can also output .mng (animated PNG, not well-supported), but ffmpeg is needed as a delegate for .mp4.\n"
    )
    extra_fp.write(
        b"#convert -verbose -delay $delay -loop 0 $output_base* -resize $size animation.gif\n"
    )
    extra_fp.write(b"#Or ffmpeg can output .mp4 directly.\n")
    extra_fp.write(b"ffmpeg -f image2 -i $output_base%d.png animation.mp4\n")
    extra_fp.close()

    # Create montage script
    extra_file = Path(plot_options.outputfolder) / "montage.sh"
    try:
        extra_fp = extra_file.open("wb")
    except OSError:
        logging.error("The montage.sh script couldn't be created.")

    extra_fp.write(b"#!/bin/bash\n")
    extra_fp.write(b"set -x\n")
    extra_fp.write(f'output_base="{plot_options.output_base}"\n'.encode())
    extra_fp.write(b"montage $output_base* -geometry +2+2 montage.png\n")
    extra_fp.close()


def run_gmt_scripts(options: Options) -> None:
    """Run GMT plotting scripts.

    On Linux, runs create_plots.sh directly via bash.

    Args:
        options: Options object. Contains:
                 - output_folder: Path to the output directory.
    """
    script_path = options.output_folder / "create_plots.sh"
    script_path.chmod(0o755)

    logging.info(f"Running {os.fspath(script_path)}...")
    result = subprocess.run(
        ["bash", script_path.name],
        cwd=os.fspath(script_path.parent),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        logging.error(f"GMT script failed (exit {result.returncode}):\n{result.stderr}")
    else:
        logging.info("GMT scripts completed successfully.")
    if result.stdout:
        logging.debug(f"GMT stdout:\n{result.stdout}")


# ---------------------------------------------------------------------------
# Extraction helpers (called from main)
# ---------------------------------------------------------------------------


def zip_script(options: Options) -> None:
    """Zip the current script into the output folder.

    Args:
        options: Options object (output_folder).
    """
    current_script_path = Path(__file__).resolve()
    zip_file_path = options.output_folder / (current_script_path.name + ".zip")
    with zipfile.ZipFile(zip_file_path, "w") as zipf:
        zipf.write(current_script_path, arcname=current_script_path.name)
    logging.info(
        f"Successfully zipped {os.fspath(current_script_path)} to {os.fspath(zip_file_path)}"
    )


def load_input_data(options: Options) -> xr.Dataset:
    """Load input data based on the configured input_choice.

    Args:
        options: Options object (input_choice, and all folder paths).

    Returns:
        Loaded xarray Dataset.
    """
    if options.use_kerchunk and options.input_choice in _KERCHUNK_CONFIG:
        return load_via_kerchunk(options)

    if options.input_choice == "SSHA":
        input_data = load_ssha_files(options)
    elif options.input_choice == "Argo":
        input_data = load_argo_file(options)
        input_data["ohc_2d_lt_700m"] = input_data["ohc_2d"] - input_data["ohc_2d_700m"]
    elif options.input_choice == "MUR_SST":
        input_data = load_mur_sst_files(options)
    elif options.input_choice == "AQUA_MODIS":
        input_data = load_aqua_modis_files(options)
        logging.info(input_data)
    elif options.input_choice == "GRACE":
        input_data = load_grace_file(options)
    elif options.input_choice == "simple_grids":
        input_data = load_simple_grid_files(options)
    else:
        logging.error(f"DID NOT RECOGNIZE {options.input_choice = }")
        raise ValueError(f"Unknown input_choice: {options.input_choice}")
    return input_data


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
