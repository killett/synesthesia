from __future__ import annotations

import argparse
import csv
import logging
import os
import shutil
import sys
import timeit
import zipfile
import datetime as dt
from pathlib import Path
from typing import Any, BinaryIO, Callable, Final

import h5py  # noqa: F401 — needed by xarray for NetCDF4 backend
import matplotlib.pyplot as plt
import nfft
import numpy as np
import pandas as pd
import statsmodels.api as sm
import xarray as xr

from colour import SpectralDistribution, sd_to_XYZ, XYZ_to_sRGB
from colour.colorimetry import MSDS_CMFS_STANDARD_OBSERVER

__version__ = "0.1.1"
DAYS_IN_YEAR: Final[float] = 365.25


# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------

class Options:
    """All global options in one place."""

    def __init__(self) -> None:
        """Initialize Options with default values."""
        # Identity
        self.my_name:           str  = Path(sys.argv[0]).stem
        self.log_mode:          int  = logging.INFO
        self.args: argparse.Namespace | None = None

        # Paths
        self.output_base:       Path = Path("./output")
        self.ssha_folder:       Path = Path("./sealevel_spectra/fast_202306/fast_netCDF4")
        self.argo_folder:       Path = Path("./sealevel_spectra/Argo")
        self.mur_sst_folder:    Path = Path("./sealevel_spectra/MUR_SST/MUR25-JPL-L4-GLOB-v04.2")
        self.aqua_modis_folder: Path = Path("./sealevel_spectra/AQUA_MODIS")
        self.grace_folder:      Path = Path("./sealevel_spectra/JPL_GRACE_mascons")
        self.cie_file:          Path = Path("./sealevel_spectra/ciexyz31_1_trimmed_400nm_700nm.csv")
        self.output_folder:     Path = Path()  # set in main()

        # Numerical knobs
        self.dpi:               int   = 300
        self.input_choice:      str   = "SSHA"
        self.xskip:             int   = 6
        self.min_period:        float = 30.0
        self.max_period:        float = 60.0
        self.thepower:          float = 0.8
        self.figsize:     tuple[int, int] = (10, 5)

        # Data key names (set per input_choice in configure_keys_for_input)
        self.x_key:             str   = "Time"
        self.y_key:             str   = "SLA"
        self.lat_key:           str   = "Latitude"
        self.lon_key:           str   = "Longitude"

        # Function selection
        self.use_new_funcs:     bool  = True

        # Files to copy into output
        self.files_to_copy: list[str] = ["projections.sh", "overflow.sh", "notation.sh"]

        # Large dicts (initialized by helpers in main)
        self.plot_options:  dict[str, Any] = {}
        self.grid:          dict[str, Any] = {}
        self.results:       dict[str, Any] = {}


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
        description=f"Spectral color mapping. {options.my_name} version {__version__}")
    parser.add_argument("-v", "--version",
                        action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("description", type=str, nargs="?", default="",
                        help="A description string encapsulated in quotes.")
    parser.add_argument("--input-choice", type=str, default=None,
                        choices=["SSHA", "Argo", "MUR_SST", "AQUA_MODIS", "GRACE"],
                        help=f"Data source (default: {options.input_choice}).")
    parser.add_argument("--xskip", type=int, default=None,
                        help="Skip every N points in lat/lon (default: per-dataset).")
    parser.add_argument("--min-period", type=float, default=None,
                        help="Minimum period in days (default: per-dataset).")
    parser.add_argument("--max-period", type=float, default=None,
                        help="Maximum period in days (default: per-dataset).")
    parser.add_argument("--dpi", type=int, default=options.dpi,
                        help=f"DPI for output images (default: {options.dpi}).")
    parser.add_argument("-d", "--debug", action="store_true",
                        help="Enable debug logging.")
    options.args = parser.parse_args()
    assert options.args is not None  # For mypy
    if options.args.debug:
        options.log_mode = logging.DEBUG
    if options.args.input_choice is not None:
        options.input_choice = options.args.input_choice
    options.dpi = options.args.dpi


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    """Main function."""
    start_time = dt.datetime.now()

    options = Options()
    parse_arguments(options)
    assert options.args is not None  # For mypy

    # Logging setup (replaces deleted logging_setup function)
    logging.basicConfig(level=options.log_mode,
                        format="%(asctime)s - %(levelname)s - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")

    # Create output folder
    date_str = start_time.strftime("%Y%m%d-%H%M%S")
    options.output_folder = options.output_base / date_str
    if options.args.description:
        options.output_folder = options.output_folder.parent / (
            options.output_folder.name + " - " + options.args.description)
    options.output_folder.mkdir(parents=True, exist_ok=True)
    if not options.output_folder.is_dir():
        raise ValueError(f"!!! Problem creating {os.fspath(options.output_folder)}")

    # Copy support files
    for file_name in options.files_to_copy:
        shutil.copy(file_name, options.output_folder)

    # Zip the current script into the output folder
    zip_script(options)

    # Initialize plot_options, grid, results dicts
    init_plot_options(options)
    init_grid(options)
    init_results(options)

    # Matplotlib setup
    plt.style.use("dark_background")
    plt.rcParams["font.size"]      = 14
    plt.rcParams["axes.linewidth"] = 2

    # Function variant selection
    if 0:
        spectrum2xyz_fn = spectrum2xyz_old
        xyz2rgb_fn      = xyz2rgb_old
        funcs_desc      = " OLD functions"
    else:
        spectrum2xyz_fn = spectrum2xyz_new
        xyz2rgb_fn      = xyz2rgb_new
        funcs_desc      = " NEW functions"

    # Configure data keys per input_choice
    configure_keys_for_input(options)

    # Backfill None CLI values with per-dataset defaults
    for attr in ("xskip", "min_period", "max_period"):
        if getattr(options.args, attr) is None:
            setattr(options.args, attr, getattr(options, attr))
        else:
            setattr(options, attr, getattr(options.args, attr))

    # Load and decimate input data
    logging.info(f"Loading {options.input_choice}")
    input_data = load_input_data(options)
    logging.info(f"Grabbing one lat/lon point in every {options.xskip ** 2} points...")
    input_data = input_data.isel({
        options.lat_key : slice(None, None, options.xskip),
        options.lon_key : slice(None, None, options.xskip),
    })
    logging.info(" done.")

    # Ensure even number of time steps for NFFT
    if input_data.sizes[options.x_key] % 2 == 1:
        logging.info(f"Deleting last data point because number of time stamps needs to be even for NFFT. Before deletion: {input_data.sizes[options.x_key] = }")
        input_data = input_data.isel({options.x_key: slice(None, -1)})

    # Validate min/max period against available NFFT periods
    sliced_data    = input_data.isel({options.lat_key: 0, options.lon_key: 0})
    power_spectrum = nfft_power(options, sliced_data)
    power_spectrum = convert_spectrum_from_frequency_to_period(power_spectrum)

    logging.info(f"{options.min_period = } and {np.min(power_spectrum.period.values) = }")
    if options.min_period < np.min(power_spectrum.period.values) or options.min_period >= np.max(power_spectrum.period.values):
        logging.error(f"!!! WARNING!!! originally {options.min_period = } but {np.min(power_spectrum.period.values) = } and {np.max(power_spectrum.period.values) = }")
        options.min_period = np.min(power_spectrum.period.values)
        logging.error(f"So now {options.min_period = } which equals {np.min(power_spectrum.period.values) = }")
    logging.info(f"{options.max_period = } and {np.max(power_spectrum.period.values) = }")
    if options.max_period <= np.min(power_spectrum.period.values) or options.max_period > np.max(power_spectrum.period.values):
        logging.error(f"!!! WARNING!!! originally {options.max_period = } but {np.min(power_spectrum.period.values) = } and {np.max(power_spectrum.period.values) = }")
        options.max_period = np.max(power_spectrum.period.values)
        logging.error(f"So now {options.max_period = } which equals {np.max(power_spectrum.period.values) = }")

    # Load CIE color matching functions
    cie = load_cie_functions(options)

    # ===================================================================
    # HOT PATH — optimized computation loop (copied verbatim from v46)
    # ===================================================================
    hot_start = timeit.default_timer()
    logging.info("Starting timeseries_to_xyz WITHOUT dask (optimized)...")

    # === PRE-COMPUTATION (outside loop) ===
    stacked       = input_data.stack(latlon=[options.lat_key, options.lon_key])
    latlon_coord  = stacked["latlon"]

    # Extract raw numpy arrays once
    times                = stacked[options.x_key].values
    data_2d              = stacked[options.y_key].values  # shape (n_times, n_points)
    n_times, n_points    = data_2d.shape

    # Time axis in days (shared across all grid points)
    x_days = (times - times[0]).astype(float) / (24 * 3600 * 1e9)

    # Ensure even length for NFFT
    N = len(x_days)
    if N % 2:
        logging.info(f"Trimming last time step for even NFFT length: {N} -> {N - 1}")
        x_days  = x_days[:-1]
        data_2d = data_2d[:-1, :]
        N      -= 1

    # Detrending design matrix: [constant, trend, accel] (shared)
    design_matrix = np.column_stack([np.ones(N), x_days, x_days**2])

    # NFFT parameters (shared)
    x_min   = x_days.min()
    x_range = x_days.max() - x_min
    x_norm  = (x_days - x_min) / x_range - 0.5
    k       = -(N // 2) + np.arange(N)
    xf      = k / x_range
    xf_half = xf[N // 2 + 1:]  # positive frequencies, ascending

    # Period array (ascending) for mapping power spectrum to wavelength
    periods_ascending  = (1.0 / xf_half)[::-1]
    cie_wavelengths    = cie["wavelength"].values
    n_wl               = len(cie_wavelengths)
    wavelength_targets = np.linspace(options.min_period, options.max_period, n_wl)

    # Pre-compute boundary handling for map_power_spectrum (replicates original reindex+nearest)
    periods_extended = periods_ascending.copy()
    if options.min_period not in periods_ascending:
        periods_extended = np.append(periods_extended, options.min_period)
    if options.max_period not in periods_ascending:
        periods_extended = np.append(periods_extended, options.max_period)
    periods_extended = np.sort(periods_extended)

    # For each extended period, find nearest original period index (for nearest-fill)
    nearest_indices = np.array([np.argmin(np.abs(periods_ascending - p)) for p in periods_extended])

    # Filter to [min_period, max_period]
    period_mask        = (periods_extended >= options.min_period) & (periods_extended <= options.max_period)
    interp_periods     = periods_extended[period_mask]
    interp_nearest_idx = nearest_indices[period_mask]

    # CIE color matching functions for sd_to_XYZ
    cmfs = MSDS_CMFS_STANDARD_OBSERVER["CIE 1931 2 Degree Standard Observer"]

    # Pre-filter all-NaN grid points (land)
    valid_mask = ~np.all(np.isnan(data_2d), axis=0)
    n_valid    = valid_mask.sum()
    logging.info(f"Processing {n_valid} valid grid points out of {n_points} total...")

    # === MAIN LOOP (pure numpy, no xarray overhead) ===
    xyz_results = np.full((n_points, 3), np.nan)

    for i in range(n_points):
        if not valid_mask[i]:
            continue

        y = data_2d[:, i]

        # 1. Detrend with numpy lstsq (replaces statsmodels OLS)
        params, _, _, _ = np.linalg.lstsq(design_matrix, y, rcond=None)
        detrended       = y - design_matrix @ params

        # 2. NFFT
        f_k   = nfft.nfft(x_norm, detrended)
        power = np.abs(f_k)**2

        # 3. Take positive frequencies, reverse to ascending period order
        power_half      = power[N // 2 + 1:]
        power_ascending = power_half[::-1]

        # 4. Map power spectrum to wavelength grid (replicates map_power_spectrum boundary handling)
        power_interp_src = power_ascending[interp_nearest_idx]
        mapped_power     = np.interp(wavelength_targets, interp_periods, power_interp_src)

        # 5. Spectrum to XYZ using colour-science sd_to_XYZ
        spd             = SpectralDistribution(mapped_power, cie_wavelengths)
        xyz_results[i]  = sd_to_XYZ(spd, cmfs)

    logging.info("Main loop complete.")

    # === POST-LOOP: Normalize XYZ ===
    max_y = np.nanmax(xyz_results[:, 1])
    logging.info(f"The highest value of 'y' is: {max_y}")
    xyz_results /= max_y

    # Raise Y to power (brightens dark areas)
    logging.info(f"Raising y to power {options.thepower} while keeping chromaticity constant. (This brightens dark areas.)")
    factor       = np.power(xyz_results[:, 1], 1.0 - options.thepower)
    xyz_results /= factor[:, np.newaxis]

    # === POST-LOOP: Vectorized XYZ -> sRGB (single call, replaces groupby.map) ===
    logging.info("Converting to RGB...")
    RGB_results = XYZ_to_sRGB(xyz_results)  # shape (n_points, 3)

    # === POST-LOOP: Vectorized fix_gamut (replaces groupby.map) ===
    logging.info("Fixing RGB out-of-gamut values and normalizing...")
    rgb_arr = RGB_results.T.copy()  # shape (3, n_points) for easier per-channel access

    nan_mask  = np.any(np.isnan(rgb_arr), axis=0)
    min_vals  = np.full(n_points, 0.0)
    min_vals[~nan_mask] = rgb_arr[:, ~nan_mask].min(axis=0)
    needs_fix = (~nan_mask) & (min_vals < 0)

    if np.any(needs_fix):
        offset     = -min_vals[needs_fix] + 1.0 / 255.0
        y_vals     = xyz_results[needs_fix, 1]
        fix_factor = y_vals / (y_vals + offset)
        rgb_arr[:, needs_fix] = (rgb_arr[:, needs_fix] + offset) * fix_factor

    # Normalize RGB
    highest_value = np.nanmax(rgb_arr)
    rgb_arr      /= highest_value

    # === Assemble xarray Dataset and unstack ===
    rgb_ds = xr.Dataset({
        "x"     : ("latlon", xyz_results[:, 0]),
        "y"     : ("latlon", xyz_results[:, 1]),
        "z"     : ("latlon", xyz_results[:, 2]),
        "red"   : ("latlon", rgb_arr[0]),
        "green" : ("latlon", rgb_arr[1]),
        "blue"  : ("latlon", rgb_arr[2]),
    }, coords={"latlon": latlon_coord})

    spectral_color_maps = xr.Dataset({key: rgb_ds[key].unstack("latlon") for key in rgb_ds.data_vars})

    # Stop the timer
    hot_end    = timeit.default_timer()
    time_taken = hot_end - hot_start
    logging.info(f"Time taken WITHOUT DASK: {time_taken:.2f} seconds")

    # Convert RGB values from [0,1] to [0,255]
    #for thiskey in ['red','green','blue']:
    #    spectral_color_maps[thiskey] = spectral_color_maps[thiskey] * 255.0

    # ===================================================================
    # Save outputs and plot
    # ===================================================================
    rgb_filenames = []
    for thekey in spectral_color_maps.keys():
        date_str = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
        filename = options.output_folder / f"{thekey.replace(' ', '_')}.nc"
        rgb_filenames.append(filename)

        logging.info(f"Saving {os.fspath(filename)}...")
        spectral_color_maps[thekey].to_netcdf(filename)
    logging.info("Finished saving.")

    for thekey in spectral_color_maps.keys():
        img = plt.imshow(spectral_color_maps[thekey], origin="lower")
        plt.colorbar(img, orientation="horizontal")
        plt.title(thekey)
        plt.savefig(options.output_folder / f"output_{thekey}.png", dpi=options.dpi)
        plt.close()

    # Stack into an RGB image
    image = np.dstack((
        spectral_color_maps["red"].values,
        spectral_color_maps["green"].values,
        spectral_color_maps["blue"].values,
    ))

    plt.imshow(image, origin="lower")
    plt.title(f"{options.y_key}, periods {options.min_period / 7.:.0f}-{options.max_period / 7.:.0f} weeks,{funcs_desc}")
    plt.savefig(options.output_folder / "image_matplotlib.png", dpi=options.dpi)
    plt.close()

    #write_gmt_scripts(options)
    #run_gmt_scripts(options)

    logging.error("!!!WARNING!!! Next line assumes these units are originally in ns and you want the units to be days!!!")
    logging.info("All finished!")
    logging.info("Download and analyze chlorophyll data, as well as sea surface salinity data!")

    end_time     = dt.datetime.now()
    elapsed_time = end_time - start_time
    logging.info(f"Elapsed time: {elapsed_time}")
    logging.shutdown()


# ---------------------------------------------------------------------------
# Pure helper functions (no options needed)
# ---------------------------------------------------------------------------

def gaussian(x: np.ndarray, mu: float, sig: float) -> np.ndarray:
    """Compute a Gaussian function."""
    return np.exp(-np.power(x - mu, 2.0) / (2 * np.power(sig, 2.0)))


def spectrum2xyz_old(spectrum: xr.Dataset, cie: xr.Dataset) -> xr.Dataset:
    """Convert a power spectrum to XYZ tristimulus values using Riemann sum integration."""
    xyz = {}
    wavelength_step_size = np.diff(cie["wavelength"].values).mean()
    for l in ["x", "y", "z"]:
        temp_values = spectrum["power"].values * cie[l].values
        xyz[l] = (temp_values.sum() * wavelength_step_size)
    return xr.Dataset(xyz)


def spectrum2xyz_new(spectrum: xr.Dataset, cie: xr.Dataset) -> xr.Dataset:
    """Convert a power spectrum to XYZ tristimulus values using colour-science.

    Args:
        spectrum: Dataset with 'power' variable and 'wavelength' coordinate.
        cie:      CIE color matching functions (unused, kept for API compatibility).

    Returns:
        Dataset with 'x', 'y', 'z' tristimulus values.
    """
    spd  = SpectralDistribution(spectrum["power"].values, spectrum.coords["wavelength"].values)
    cmfs = MSDS_CMFS_STANDARD_OBSERVER["CIE 1931 2 Degree Standard Observer"]
    XYZ  = sd_to_XYZ(spd, cmfs)
    return xr.Dataset({"x": XYZ[0], "y": XYZ[1], "z": XYZ[2]})


def raise_y_to_power(xyz: xr.Dataset, power: float) -> xr.Dataset:
    """Raise Y to a power while keeping chromaticity constant."""
    power  = 1 - power
    factor = pow(xyz["y"].values, power)
    for key in xyz.data_vars:
        xyz[key].values /= factor
    return xyz


def xyz2rgb_old(xyz: xr.Dataset) -> xr.Dataset:
    """Convert XYZ to sRGB using a manual matrix multiply."""
    A = np.array([[3.2409699, -1.5373832, -0.49861079],
                  [-0.96924375, 1.8759676, 0.041555082],
                  [0.055630032, -0.20397685, 1.0569714]])

    xyz_vector = np.array([xyz["x"], xyz["y"], xyz["z"]])
    rgb        = np.dot(A, xyz_vector)
    rgb        = {"red": rgb[0], "green": rgb[1], "blue": rgb[2]}
    result     = xr.merge([xyz, rgb])
    return result


def xyz2rgb_new(xyz: xr.Dataset) -> xr.Dataset:
    """Convert XYZ to sRGB using colour-science."""
    XYZ    = np.array([xyz["x"].values.squeeze(), xyz["y"].values.squeeze(), xyz["z"].values.squeeze()])
    RGB    = XYZ_to_sRGB(XYZ)
    rgb    = xr.Dataset({"red": RGB[0], "green": RGB[1], "blue": RGB[2]})
    result = xr.merge([xyz, rgb])
    return result


def fix_gamut(rgb: xr.Dataset) -> xr.Dataset:
    """Fix out-of-gamut RGB values by adding white light to make all values positive."""
    R          = rgb["red"].values.squeeze()
    G          = rgb["green"].values.squeeze()
    B          = rgb["blue"].values.squeeze()
    rgb_values = np.array([R, G, B])

    if np.any(np.isnan(rgb_values)):
        return rgb

    min_val = rgb_values.min()
    if min_val < 0:
        # Make "min" positive as in paper's appendix, and add 1/255 so the rescale value isn't 0.
        min_val    = -min_val + 1.0 / 255.0
        # This factor rescales the luminance back to its original value.
        factor     = rgb["y"].values.squeeze() / (rgb["y"].values.squeeze() + min_val)
        rgb_values = (rgb_values + min_val) * factor

    result = xr.Dataset({
        "x"     : rgb["x"],
        "y"     : rgb["y"],
        "z"     : rgb["z"],
        "red"   : rgb_values[0],
        "green" : rgb_values[1],
        "blue"  : rgb_values[2],
    })
    return result


def gamma_correct_rgb(rgb: xr.Dataset) -> xr.Dataset:
    """Apply gamma correction to RGB values."""
    gamma_inv = 0.45
    crit      = 0.018  # RGB values are gamma corrected differently below and above crit.
    h         = 4.506813168
    g         = -0.09914989
    f         = 1.09914989
    for key in ["red", "green", "blue"]:
        # Typo in Hughes and Williams 2010 equation A7, compared to Charles Poynton's GammaFAQ:
        # http://www.poynton.com/GammaFAQ.html
        if rgb[key] <= crit:
            rgb[key] *= h
        else:
            rgb[key] = f * pow(rgb[key], gamma_inv) + g
    return rgb


def rms_and_mean(x: xr.DataArray) -> xr.Dataset:
    """Compute RMS and mean of an array."""
    rms_value  = np.sqrt(np.mean(x**2))
    mean_value = np.mean(x)
    return xr.Dataset({"rms": rms_value, "mean": mean_value})


def normalize_data(data: Any, max_value: float) -> list[float]:
    """Normalize data by dividing by max_value."""
    return [i / max_value for i in data]


# ---------------------------------------------------------------------------
# GMT pure helpers (no options needed)
# ---------------------------------------------------------------------------

def write_gmt_defs(new_fp: BinaryIO) -> None:
    """Write clarifying GMT definitions to the script file."""
    new_fp.write(b"#######################################################\n")
    new_fp.write(b"#Clarifying definitions. Do not change!################\n")
    new_fp.write(b"start=\" -K \" #Should always redirect using > to write new PS.\n")
    new_fp.write(b"middle=\" -O -K \" #Should always redirect using >> to append to PS.\n")
    new_fp.write(b"end=\" -O \" #Should always redirect using >> to append to PS.\n")
    new_fp.write(b"#######################################################\n")


def write_gmt_coastlines(new_fp: BinaryIO) -> None:
    """Write GMT coastlines commands to the script file."""
    new_fp.write(b"gmt coast -W$coast_thk/$coast_color $coast_res $range $projection $map_pos $middle >> $plot_base.ps\n")
    new_fp.write(b"if [ $coastlines == 2 ]\n")
    new_fp.write(b"then\n")
    new_fp.write(b"  gmt plot -N $coast_file -: -Sc$coast_thk -W$coast_thk/$coast_color $range $projection $map_pos $middle >> $plot_base.ps\n")
    new_fp.write(b"fi\n")


def finish_flip_backgrounds(flip_fp: BinaryIO) -> None:
    """Finish writing the flip_backgrounds.sh script."""
    flip_fp.write(b"for current_base in $all_bases\n")
    flip_fp.write(b"do\n")
    flip_fp.write(b"  #Change black to cyan temporarily.\n")
    flip_fp.write(b"  convert $current_base.png -fill cyan -opaque black $current_base.png\n")
    flip_fp.write(b"  #Change white to black.\n")
    flip_fp.write(b"  convert $current_base.png -fill black -opaque white $current_base.png\n")
    flip_fp.write(b"  #Change temporary cyan to white.\n")
    flip_fp.write(b"  convert $current_base.png -fill white -opaque cyan $current_base.png\n")
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
    wavelengths: list[float]     = []

    try:
        with cie_path.open("r") as in_fp:
            reader = csv.reader(in_fp)
            for row in reader:
                wavelengths.append(float(row[0]))
                for i, val in enumerate(row[1:], start=0):
                    data[list(data.keys())[i]].append(float(val))
    except IOError:
        logging.error(f"The CIE file, {os.fspath(cie_path)}, failed to open.")

    da = xr.Dataset(
        {var: ("wavelength", data[var]) for var in data},
        coords={"wavelength": wavelengths},
    )
    da.attrs["title"]                     = "CIE 1931 color matching functions"
    da.coords["wavelength"].attrs["units"] = "nm"
    return da


def synthetic_spectrum(cie: xr.Dataset, mu: float, sig: float) -> xr.Dataset:
    """Create a synthetic Gaussian spectrum on the CIE wavelength grid."""
    wavelengths = cie.coords["wavelength"].values
    power       = gaussian(wavelengths, mu, sig)
    spectrum    = xr.DataArray(power, coords=[("wavelength", wavelengths)], name="power").to_dataset()
    return spectrum


def synthetic_timeseries(options: Options, signal: str = "annual", signal_amplitude: float = 1,
                         noise: str = "white", noise_level: float = 0.1,
                         temporal_resolution: str = "monthly",
                         time_start: dt.datetime = dt.datetime(2001, 1, 1),
                         time_stop: dt.datetime = dt.datetime(2005, 1, 1)) -> xr.Dataset:
    """Generate a synthetic timeseries with signal and noise.

    Args:
        options:              Options object. Contains:
                                  - x_key: Time coordinate name.
                                  - y_key: Variable name.
        signal:               Signal type ("annual").
        signal_amplitude:     Amplitude of the signal.
        noise:                Noise type ("white").
        noise_level:          Standard deviation of noise.
        temporal_resolution:  "monthly" or "daily".
        time_start:           Start date.
        time_stop:            End date.

    Returns:
        Dataset with the synthetic timeseries.
    """
    if temporal_resolution == "monthly":
        dates = pd.date_range(start=time_start, end=time_stop, freq="M") + pd.Timedelta(days=15)
    elif temporal_resolution == "daily":
        dates = pd.date_range(start=time_start, end=time_stop, freq="D")
    else:
        logging.error(f"!!!WARNING!!! {temporal_resolution = }")

    if signal == "annual":
        t             = (dates - time_start).days / 365.25
        signal_values = signal_amplitude * np.sin(2 * np.pi * t)
    else:
        logging.error(f"!!!WARNING!!! {signal = }")

    if noise == "white":
        noise_values = noise_level * np.random.randn(len(dates))
    else:
        logging.error(f"!!!WARNING!!! {noise = }")

    measurements = signal_values + noise_values

    return xr.Dataset(
        {options.y_key: (options.x_key, measurements)},
        coords={options.x_key: dates},
    )


def fancy_detrend(timeseries: xr.Dataset, x_key: str, y_key: str,
                  terms: list[str] | None = None) -> tuple[xr.Dataset, dict[str, float]]:
    """Detrend a timeseries using OLS regression.

    Args:
        timeseries: Input timeseries dataset.
        x_key:      Name of the time coordinate.
        y_key:      Name of the data variable.
        terms:      Regression terms (e.g. ["constant", "trend", "accel"]).

    Returns:
        Tuple of (detrended timeseries, dict of fit coefficients).
    """
    if terms is None:
        terms = ["constant", "trend"]

    x = (timeseries[x_key] - timeseries[x_key][0]).values.astype(float) / (24 * 3600 * 1e9)
    y = timeseries[y_key].values.squeeze()

    design_matrix = []

    if "constant" in terms:
        design_matrix.append(np.ones_like(x))
    if "trend" in terms:
        design_matrix.append(x)
    if "accel" in terms:
        design_matrix.append(x**2)
    if "annual" in terms:
        design_matrix.append(np.sin(2 * np.pi * x / DAYS_IN_YEAR))
        design_matrix.append(np.cos(2 * np.pi * x / DAYS_IN_YEAR))
    if "semiannual" in terms:
        design_matrix.append(np.sin(2 * np.pi * x / DAYS_IN_YEAR / 2.0))
        design_matrix.append(np.cos(2 * np.pi * x / DAYS_IN_YEAR / 2.0))

    design_matrix = np.column_stack(design_matrix)

    model  = sm.OLS(y, design_matrix)
    result = model.fit()

    detrended_y = y - result.fittedvalues

    detrended_timeseries             = timeseries.copy()
    detrended_timeseries[y_key]      = (timeseries[y_key].dims, detrended_y.reshape(timeseries[y_key].shape))

    fits: dict[str, float] = {}
    for i, term in enumerate(terms):
        if term == "annual":
            fits["annual_sin"] = result.params[i]
            fits["annual_cos"] = result.params[i + 1]
        elif term == "semiannual":
            fits["semiannual_sin"] = result.params[i]
            fits["semiannual_cos"] = result.params[i + 1]
        else:
            fits[term] = result.params[i]

    return detrended_timeseries, fits


def turn_fits_into_timeseries(timeseries: xr.Dataset, x_key: str, y_key: str,
                              fits: dict[str, float]) -> xr.Dataset:
    """Reconstruct a timeseries from fit coefficients.

    Args:
        timeseries: Original timeseries (provides time axis and shape).
        x_key:      Name of the time coordinate.
        y_key:      Name of the data variable.
        fits:       Dict of fit coefficients from fancy_detrend.

    Returns:
        Dataset with reconstructed fitted values.
    """
    x             = (timeseries[x_key] - timeseries[x_key][0]).values.astype(float) / (24 * 3600 * 1e9)
    fitted_values = np.zeros_like(x)

    for term, fit_value in fits.items():
        if term == "constant":
            fitted_values += fit_value
        elif term == "trend":
            fitted_values += fit_value * x
        elif term == "accel":
            fitted_values += fit_value * x**2
        elif term == "annual_sin":
            fitted_values += fit_value * np.sin(2 * np.pi * x / DAYS_IN_YEAR)
        elif term == "annual_cos":
            fitted_values += fit_value * np.cos(2 * np.pi * x / DAYS_IN_YEAR)
        elif term == "semiannual_sin":
            fitted_values += fit_value * np.sin(2 * np.pi * x / DAYS_IN_YEAR / 2.0)
        elif term == "semiannual_cos":
            fitted_values += fit_value * np.cos(2 * np.pi * x / DAYS_IN_YEAR / 2.0)

    fitted_timeseries        = timeseries.copy()
    fitted_timeseries[y_key] = (timeseries[y_key].dims, fitted_values.reshape(timeseries[y_key].shape))
    return fitted_timeseries


def nfft_power(options: Options, timeseries: xr.Dataset) -> xr.Dataset:
    """Compute the power spectrum of a timeseries using NFFT.

    Args:
        options:    Options object. Contains:
                        - x_key: Time coordinate name.
                        - y_key: Variable name.
        timeseries: Input timeseries dataset.

    Returns:
        Dataset with 'power' variable and 'frequency' coordinate.
    """
    #logging.error("!!!WARNING!!! Next line assumes these units are originally in ns and you want the units to be days!!!")
    x = (timeseries[options.x_key] - timeseries[options.x_key][0]).values.astype(float).squeeze() / (24 * 3600 * 1e9)
    y = timeseries[options.y_key].values.squeeze()

    # If timeseries has an odd number of points,
    # remove the last data point, then calculate min, range.
    N = 1
    while N % 2:
        N       = len(x)
        x_min   = np.min(x)
        x_range = np.max(x) - np.min(x)
        x_norm  = (x - x_min) / x_range - 0.5
        if N % 2:
            logging.error(f"!!! WARNING!!! LENGTH NEEDS TO BE EVEN FOR NFFT, BUT: {len(x) = }")
            logging.error("!!! DELETING LAST DATA POINT!")
            x = np.delete(x, -1)
            y = np.delete(y, -1)

    # Define Fourier modes
    k  = -(N // 2) + np.arange(N)
    # Convert Fourier modes to frequencies
    xf = k / x_range

    # Perform NFFT.
    f_k = nfft.nfft(x_norm, y)

    # Compute power spectrum, which is the square of the absolute value of the Fourier Transform
    power_spectrum = np.abs(f_k)**2

    # Only take the positive frequencies. Since the output is symmetric, this will not lose any information.
    power_spectrum = power_spectrum[N // 2 + 1:]
    xf_half        = xf[N // 2 + 1:]

    # Create xarray DataArray with coordinates
    power_spectrum_da = xr.DataArray(power_spectrum, coords=[("frequency", xf_half)], name="power")

    # Convert this DataArray to a Dataset
    spectrum                              = power_spectrum_da.to_dataset()
    spectrum["frequency"].attrs["units"]  = "1/days"

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
    new_spectrum = new_spectrum.swap_dims({"frequency": "period"}).drop_vars("frequency")
    new_spectrum = new_spectrum.sortby("period")

    return new_spectrum


def map_power_spectrum(cie: xr.Dataset, power_spectrum: xr.Dataset,
                       min_period: float = -1, max_period: float = -1) -> xr.Dataset:
    """Map a power spectrum from period space to the CIE wavelength grid.

    Args:
        cie:             CIE color matching functions dataset.
        power_spectrum:  Power spectrum with 'period' coordinate.
        min_period:      Minimum period for mapping.
        max_period:      Maximum period for mapping.

    Returns:
        Dataset with 'power' mapped onto the CIE 'wavelength' coordinate.
    """
    existing_periods = power_spectrum["period"].values

    new_period_added = False

    if min_period not in existing_periods:
        existing_periods = np.append(existing_periods, min_period)
        new_period_added = True
    if max_period not in existing_periods:
        existing_periods = np.append(existing_periods, max_period)
        new_period_added = True

    if new_period_added:
        new_periods    = np.sort(existing_periods)
        power_spectrum = power_spectrum.reindex(period=new_periods, method="nearest")

    power_spectrum["power"] = power_spectrum["power"].interpolate_na(dim="period")

    power_spectrum = power_spectrum.where(
        (power_spectrum["period"] >= min_period) & (power_spectrum["period"] <= max_period), drop=True)

    mapped_power_values = np.interp(
        np.linspace(min_period, max_period, len(cie["wavelength"])),
        power_spectrum["period"],
        power_spectrum["power"],
    )

    new_power_spectrum = xr.Dataset(
        {"power": (("wavelength",), mapped_power_values)},
        coords={"wavelength": cie["wavelength"]},
    )
    return new_power_spectrum


# ---------------------------------------------------------------------------
# Plot functions
# ---------------------------------------------------------------------------

def plot_timeseries(options: Options, ds: xr.Dataset, title: str) -> None:
    """Plot a timeseries and save to the output folder.

    Args:
        options: Options object. Contains:
                     - figsize:       Figure size.
                     - x_key:         Time coordinate name.
                     - y_key:         Variable name.
                     - output_folder: Output directory.
                     - dpi:           DPI for saved image.
        ds:      Timeseries dataset.
        title:   Plot title.
    """
    plt.figure(figsize=options.figsize)
    plt.plot(ds[options.x_key], ds[options.y_key], color="lime")
    plt.scatter(ds[options.x_key], ds[options.y_key], marker="s", color="cyan", s=10)
    plt.title(title, color="white")
    date_str = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    filename = options.output_folder / f"{date_str}_{title.replace(' ', '_')}.png"
    plt.savefig(filename, dpi=options.dpi, format="png", transparent=False, bbox_inches="tight", facecolor="black")


def plot_fft_spectrum(options: Options, power_spectrum: xr.Dataset, title: str) -> None:
    """Plot an FFT power spectrum and save to the output folder.

    Args:
        options:        Options object (figsize, output_folder, dpi).
        power_spectrum: Power spectrum dataset with 'period' and 'power'.
        title:          Plot title.
    """
    plt.figure(figsize=options.figsize)
    plt.plot(power_spectrum.period, power_spectrum.power, color="lime", linewidth=2.0)
    plt.scatter(power_spectrum.period, power_spectrum.power, marker="s", color="cyan", s=10)
    plt.title(title)
    plt.xlabel("Period (days)")
    plt.ylabel("Power")
    date_str = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    filename = options.output_folder / f"{date_str}_{title.replace(' ', '_')}.png"
    plt.savefig(filename, dpi=options.dpi, format="png", transparent=False, bbox_inches="tight", facecolor="black")


def plot_light_spectrum(options: Options, power_spectrum: xr.Dataset, title: str) -> None:
    """Plot a light spectrum and save to the output folder.

    Args:
        options:        Options object (figsize, output_folder, dpi).
        power_spectrum: Spectrum dataset with 'wavelength' and 'power'.
        title:          Plot title.
    """
    plt.figure(figsize=options.figsize)
    plt.plot(power_spectrum.wavelength, power_spectrum.power, color="lime", linewidth=2.0)
    plt.scatter(power_spectrum.wavelength, power_spectrum.power, marker="s", color="cyan", s=10)
    plt.title(title)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Power")
    date_str = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    filename = options.output_folder / f"{date_str}_{title.replace(' ', '_')}.png"
    plt.savefig(filename, dpi=options.dpi, format="png", transparent=False, bbox_inches="tight", facecolor="black")


def plot_color(options: Options, rgb: xr.Dataset, filename: str | os.PathLike[str]) -> None:
    """Plot a solid color square and save to a file.

    Args:
        options:  Options object (dpi).
        rgb:      Dataset with 'red', 'green', 'blue' values.
        filename: Output file path.
    """
    filename = Path(filename)
    fig, ax  = plt.subplots(1, 1, figsize=(2, 2), dpi=options.dpi)
    ax.set_facecolor(tuple(rgb[key] for key in ["red", "green", "blue"]))
    ax.axis("off")
    plt.savefig(filename, dpi=options.dpi, format="png", transparent=False, bbox_inches="tight")


# ---------------------------------------------------------------------------
# Data loading functions
# ---------------------------------------------------------------------------

def load_ssha_files(options: Options, tskip: int = 1) -> xr.Dataset:
    """Load SSHA NetCDF files.

    Args:
        options: Options object. Contains:
                     - ssha_folder: Path to SSHA data directory.
        tskip:   Load every tskip-th file.

    Returns:
        Combined dataset.
    """
    sshafiles  = sorted(options.ssha_folder.glob("*.nc"))
    tskip_files = sshafiles[::tskip]
    logging.info(f"Loading {len(tskip_files)} SSHA files...")
    input_data = xr.open_mfdataset(tskip_files, data_vars="all", combine="by_coords")
    return input_data


def load_argo_file(options: Options) -> xr.Dataset:
    """Load a single Argo NetCDF file.

    Args:
        options: Options object (argo_folder).

    Returns:
        Argo dataset.
    """
    thefiles = sorted(options.argo_folder.glob("*.nc"))
    thefile  = thefiles[0]
    return xr.open_dataset(thefile)


def load_grace_file(options: Options) -> xr.Dataset:
    """Load a single GRACE NetCDF file.

    Args:
        options: Options object (grace_folder).

    Returns:
        GRACE dataset.
    """
    thefiles = sorted(options.grace_folder.glob("*.nc"))
    thefile  = thefiles[0]
    return xr.open_dataset(thefile)


def load_mur_sst_files(options: Options, tskip: int = 1) -> xr.Dataset:
    """Load MUR SST NetCDF files.

    Args:
        options: Options object (mur_sst_folder).
        tskip:   Load every tskip-th file.

    Returns:
        Combined dataset.
    """
    thefiles    = sorted(options.mur_sst_folder.glob("*2003*.nc"))
    tskip_files = thefiles[::tskip]
    logging.info(f"Loading {len(tskip_files)} MUR SST files...")
    input_data  = xr.open_mfdataset(tskip_files, data_vars="all", combine="by_coords")
    return input_data


def load_aqua_modis_files(options: Options, tskip: int = 1) -> xr.Dataset:
    """Load AQUA MODIS NetCDF files.

    Args:
        options: Options object (aqua_modis_folder).
        tskip:   Load every tskip-th file.

    Returns:
        Combined dataset.
    """
    thefiles    = sorted(options.aqua_modis_folder.glob("*2003*.nc"))
    tskip_files = thefiles[::tskip]
    logging.info(f"Loading {len(tskip_files)} AQUA_MODIS files...")
    input_data  = xr.open_mfdataset(tskip_files, data_vars="all", combine="by_coords")
    return input_data


# ---------------------------------------------------------------------------
# Analysis functions
# ---------------------------------------------------------------------------

def extract_ssha_timeseries(options: Options, ds: xr.Dataset,
                            lat: float = 30, lon: float = 135) -> xr.Dataset:
    """Extract a single-point SSHA timeseries at given coordinates.

    Args:
        options: Options object (x_key, y_key).
        ds:      Full SSHA dataset.
        lat:     Latitude (degrees).
        lon:     Longitude (degrees, converted to 0-360 if negative).

    Returns:
        Single-point timeseries dataset.
    """
    logging.info(f"Extracting SSHA timeseries at {lat = } and {lon = }")
    if lon < 0:
        lon = 360 + lon

    measurements = ds[options.y_key].sel(Latitude=lat, Longitude=lon, method="nearest").values

    ds = xr.Dataset(
        {options.y_key: (options.x_key, measurements)},
        coords={options.x_key: ds[options.x_key].values},
    )
    return ds


def timeseries_to_xyz(options: Options, timeseries: xr.Dataset, x_key: str, y_key: str,
                       min_period: float, max_period: float, cie: xr.Dataset,
                       spectrum2xyz_fn: Callable[..., xr.Dataset] = spectrum2xyz_new) -> xr.Dataset:
    """Convert a timeseries to XYZ tristimulus values via spectral analysis.

    Args:
        options:          Options object (lat_key, lon_key).
        timeseries:       Input timeseries dataset.
        x_key:            Time coordinate name.
        y_key:            Data variable name.
        min_period:       Minimum period for spectrum mapping.
        max_period:       Maximum period for spectrum mapping.
        cie:              CIE color matching functions.
        spectrum2xyz_fn:  Function to convert spectrum to XYZ.

    Returns:
        Dataset with 'x', 'y', 'z' tristimulus values.
    """
    if 0:
        lat = timeseries[options.lat_key].values
        lon = timeseries[options.lon_key].values
        #spectrum = synthetic_spectrum(cie, 500+lat, 5)
        spectrum = synthetic_spectrum(cie, 580 + 2 * lat, 5)
        #spectrum = synthetic_spectrum(cie, 550, 5)
        xyz = spectrum2xyz_fn(spectrum, cie)
        xyz["y"] = lon * lon * lon * lon
        return xyz

    fit_terms  = []
    fit_terms += ["constant", "trend", "accel"]
    #fit_terms += ['annual']
    #fit_terms += ['semiannual']
    #fit_terms += ['annual','semiannual']
    timeseries, fits = fancy_detrend(timeseries, x_key, y_key, terms=fit_terms)

    # Perform non-uniform FFT to get power spectrum.
    power_spectrum = nfft_power(options, timeseries)
    power_spectrum = convert_spectrum_from_frequency_to_period(power_spectrum)

    mapped_spectrum = map_power_spectrum(cie, power_spectrum, min_period=min_period, max_period=max_period)

    return spectrum2xyz_fn(mapped_spectrum, cie)


# ---------------------------------------------------------------------------
# GMT high-level functions
# ---------------------------------------------------------------------------

def is_polar(options: Options) -> int:
    """Determine if the data is global, north-polar, or south-polar.

    Args:
        options: Options object. Contains:
                     - results: Results dict (mutated with min/maxlat/lon).
                     - grid:    Grid dict.

    Returns:
        0 for global, 1 for north pole, 2 for south pole.
    """
    results = options.results
    grid    = options.grid
    polar   = 0

    if results["options"]["output_choice"] == 1 or results["options"]["output_choice"] == 4:
        results["minlat"] = min(grid["lat"])
        results["maxlat"] = max(grid["lat"])
        results["minlon"] = min(grid["lon"])
        results["maxlon"] = max(grid["lon"])
    elif results["options"]["output_choice"] == 5:
        results["minlat"] = min(results["latlon"]["lat"])
        results["maxlat"] = max(results["latlon"]["lat"])
        results["minlon"] = min(results["latlon"]["lon"])
        results["maxlon"] = max(results["latlon"]["lon"])
    else:
        logging.error(f"!!!!WARNING!!!!!! results['options']['output_choice'] {results['options']['output_choice']} isn't recognized.")

    if results["minlat"] < 0 and results["maxlat"] > 0:
        polar = 0
    elif results["minlat"] > 0 and results["maxlat"] > 0:
        polar = 1
    elif results["minlat"] < 0 and results["maxlat"] < 0:
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
    rgb     = {}  # Just to init choice and maxes.

    if len(results["rgb"]) != 3:
        new_fp.write(b"cpt_name=\"-Cmap.cpt \"\n")
    else:
        new_fp.write(b"cpt_name=\"-C../../rgb00001.cpt \"\n")

    if len(results["rgb"]) != 3 or results["rgb_choice"] < 2:
        if kml_output:
            new_fp.write(b"gmt colorbar $cpt_name -L $scale_format $overflow $scale_pos -A $start > scale_$plot_base.ps\n")
            new_fp.write(b"# Print units manually, otherwise they're too close to numbers on scale.\n")
            new_fp.write(b"echo $units_format $scale_units | gmt pstext -N $units_pos $misc_range $end >> scale_$plot_base.ps\n")
        else:
            new_fp.write(b"gmt colorbar $cpt_name -L $scale_format $overflow $scale_pos -A $middle >> $plot_base.ps\n")
            new_fp.write(b"# Print units manually, otherwise they're too close to numbers on scale.\n")
            new_fp.write(b"echo $units_format $scale_units | gmt pstext -N $units_pos $misc_range $middle >> $plot_base.ps\n")
    else:
        new_fp.write(f"numwidths={results['max_widths']}\n".encode())
        new_fp.write(b"gmt set TICK_LENGTH 0.3c\n")
        new_fp.write(b"scale_width=$(bc <<< \"scale=5; $scale_width / $numwidths\")\n")
        new_fp.write(b"for (( j=1; j <= $numwidths; j++ ))\n")
        new_fp.write(b"do\n")
        new_fp.write(b"  j_string=$(printf '%%05d' $j)\n")
        new_fp.write(b"  cpt_name=\"-C../rgb$j_string.cpt\"\n")
        new_fp.write(b"  gmt colorbar $cpt_name -L $scale_format $overflow $scale_pos -S -A $middle >> $plot_base.ps\n")
        new_fp.write(b"  # Print units manually, otherwise they're too close to numbers on scale.\n")
        new_fp.write(b"  echo $units_format $scale_units | gmt pstext -N $units_pos $misc_range $middle >> $plot_base.ps\n")
        new_fp.write(b"  # Move next scale to the right and get rid of the tick marks.\n")
        new_fp.write(b"  scale_x=$(bc <<< \"scale=5; $scale_x+$scale_width\")\n")
        new_fp.write(b"  scale_pos=\" -D${scale_x}c/${scale_y}c/${scale_length}c/${scale_width}c \"\n")
        new_fp.write(b"  gmt set TICK_LENGTH 0.0\n")
        new_fp.write(b"done\n")


def write_gmt_map_data(options: Options, new_fp: BinaryIO, title: str, i: int) -> None:
    """Write GMT data plotting commands to the script file.

    Args:
        options: Options object (results, grid, plot_options).
        new_fp:  File handle to the GMT script.
        title:   Map title string. If blank, outputs KMZ format.
        i:       Index of the parameter being mapped.
    """
    results      = options.results
    grid         = options.grid
    plot_options = options.plot_options

    polar      = is_polar(options)
    coastlines = plot_options["coastlines"]  # Default; overridden below for polar == 1

    if i == 0:
        new_fp.write(b"#Set resolution, coast_file, coast_thickness, and coastlines\n")
        new_fp.write(b"#on first map only because they should be universal.\n")
        coast_file = plot_options["outputfolder"] + "data/ancillary/Rignot/InSAR_GL_Antarctica.txt"
        new_fp.write(f'coast_file="{coast_file}"\n'.encode())

        if results["options"]["output_choice"] == 5:
            delta_lat = abs(results["latlon"]["lat"][1] - results["latlon"]["lat"][0])
            if delta_lat < 0.4:
                new_fp.write(b'resolution=" -E50 " #50/2000 is low/high quality.\n')
                new_fp.write(b'coast_res=" -Df+ "\n')
            elif delta_lat < 0.9:
                new_fp.write(b'resolution=" -E50 " #50/2000 is low/high quality.\n')
                new_fp.write(b'coast_res=" -Df+ "\n')
            else:
                new_fp.write(b'resolution=" -E50 " #50/2000 is low/high quality.\n')
                new_fp.write(b'coast_res=" -Di+ "\n')
        elif results["options"]["output_choice"] in [1, 4]:
            new_fp.write(b'resolution=" -E50 " #50/2000 is low/high quality.\n')
            new_fp.write(b'coast_res=" -Di+ "\n')
        else:
            logging.error(f"!!!!WARNING!!!!!! results['options']['output_choice'] {results['options']['output_choice']} isn't recognized.")

        new_fp.write(b'coast_res_orig=$coast_res #Don\'t want USA maps to repeatedly add -N2.\n')
        new_fp.write(b'coast_thk="0.6"\n')
        new_fp.write(b'coast_thk="0.009"\n')
        if polar == 1:
            coastlines = 1  # InSAR is only in Antarctica, so disable for NP plots.
        new_fp.write(f'coastlines={coastlines} #1:coast, 2:coast+InSAR.\n'.encode())

    new_fp.write(b"#coast_color is gray82 for off-white, or gray10 for dark coastlines.\n")

    # Adjust max/min latitudes for mapping points
    if results["options"]["output_choice"] in [1, 4]:
        buffer = 5
        if results["maxlat"] <= 90 - buffer:
            results["maxlat"] += buffer
        if results["minlat"] >= -90 + buffer:
            results["minlat"] -= buffer

    new_fp.write(b'title_format="0 0 30 0 0 MC"\n')
    new_fp.write(b'blurb_format="0 0 15 0 1 ML"\n')
    new_fp.write(b'units_format="0 0 13 0 0 MC"\n')
    new_fp.write(f"scale_units=\"{results['units'][i]}\"\n".encode())

    new_fp.write(b'misc_range=" -R0/1/0/1 -JX1c "\n')
    new_fp.write(b"#grdcut requires actual limits, but if grdimage uses them: GMT Fatal Error: grdimage could not allocate memory [21.69 Gb, n_items = 5823567396]\n")
    new_fp.write(b'minlon=%.3f\n' % 0.0)
    new_fp.write(b'maxlon=%.3f\n' % 360.0)
    new_fp.write(b'minlat=%.3f\n' % -90.0)
    new_fp.write(b'maxlat=%.3f\n' % 90.0)

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
            new_fp.write(f"{'#' if i > 0 else ''}projection_choice={plot_options['projection']}\n".encode())
        else:
            if polar == 1:
                new_fp.write(f"{'#' if i > 0 else ''}projection_choice=101\n".encode())
            else:
                new_fp.write(f"{'#' if i > 0 else ''}projection_choice=102\n".encode())

        new_fp.write(b"standard_circle=0 #1=all specific regions use standard circular projection.\n")
        new_fp.write(b"standard_rect=0 #1=all specific regions use standard rectangular projection.\n")
        new_fp.write(b". ./projections.sh\n")
        new_fp.write(b"if [ $projection_choice == 101 ]\n")
        new_fp.write(b"then\n")

        minlat = results["minlat"] if polar == 1 else 0.0
        new_fp.write(f"  minlat={minlat:.3f}\n".encode())
        new_fp.write(b"  actual_range=\" -R0.0/360.0/$minlat/90.0 \"\n")
        polar_radius = 90 - minlat
        new_fp.write(f"  polar_radius={polar_radius}\n".encode())
        new_fp.write(b"  projection=\" -JE0/90.0/${polar_radius}/${map_width}c \" #N. Azimuthal Equidistant\n")
        new_fp.write(b"elif [ $projection_choice == 102 ]\n")
        new_fp.write(b"then\n")

        maxlat = results["maxlat"] if polar == 2 else 0.0
        new_fp.write(f"  maxlat={maxlat:.3f}\n".encode())
        new_fp.write(b"  actual_range=\" -R0.0/360.0/-90.0/$maxlat \"\n")
        polar_radius = 90 + maxlat
        new_fp.write(f"  polar_radius={polar_radius}\n".encode())
        new_fp.write(b"  projection=\" -JE0/-90.0/{polar_radius}/${map_width}c \" #S. Azimuthal Equidistant\n")
        new_fp.write(b"fi\n")

        new_fp.write(b"range=\" -R${minlon}/${maxlon}/${minlat}/${maxlat} \"\n")
        new_fp.write(b"map_pos=\" -Xa${map_x}c -Ya${map_y}c \"\n")

        if len(results["rgb"]) == 3 and results["rgb_choice"] >= 2:
            new_fp.write(b"scale_width=1.2 #Override for RGB maps.\n")

        new_fp.write(b"scale_pos=\" -D${scale_x}c/${scale_y}c/${scale_length}c/${scale_width}c \"\n")
        new_fp.write(b"units_x=$(bc <<< \"scale=5; $scale_x+$scale_width/2\")\n")
        new_fp.write(b"units_y=$(bc <<< \"scale=5; $scale_y+$scale_length/2\")\n")
        new_fp.write(b"units_pos=\" -Xa${units_x}c -Ya${units_y}c \"\n")
        new_fp.write(b"blurb_pos=\" -Xa${blurb_x}c -Ya${blurbs_y}c \"\n")
        new_fp.write(b"blurb2_pos=\" -Xa${blurb2_x}c -Ya${blurbs_y}c \"\n")
    else:
        logging.error("!!!WARNING!!! NO TITLE!")

    # Plot data, with title on top.
    new_fp.write(f"title=\"{title}\"\n".encode())
    if len(results["rgb"]) == 3 and len(results["latlon"]["outputs"]) == 3:
        new_fp.write(b"gmt grdimage red.nc green.nc blue.nc $boundary $resolution $range $projection $map_pos $start > $plot_base.ps\n")
    else:
        new_fp.write(b"gmt grdimage $data_name $boundary $resolution $range $projection $map_pos -Cmap.cpt $start > $plot_base.ps\n")

    write_gmt_coastlines(new_fp)

    if (len(results["marker_lats"]) == len(results["latlon"]["outputs"])
            and len(results["marker_lons"]) == len(results["latlon"]["outputs"])
            and results["latlon"]["outputs"]):
        if results["marker_lats"][i] and results["marker_lons"][i]:
            new_fp.write(b"gmt plot -N $data_name -bcmarker_lons/marker_lats -S+0.5c -W5/244/164/96 -G244/164/96 $range $projection $map_pos $middle >> $plot_base.ps\n")

    new_fp.write(b"#Uncomment to put a marker at echoed coords, given as lon lat:\n")
    new_fp.write(b"#echo -85.19 -77.36 | gmt plot -N -S+0.5c -W5/244/164/96 -G244/164/96 $range $projection $map_pos $middle >> $plot_base.ps\n")

    if plot_options["plot_mascons"] != 0 and results["latlon"]["mascon_lats"]:
        new_fp.write(b"gmt plot $data_name -bcmascon_lons/mascon_lats -Sc0.01c -G139/69/19 $range $projection $map_pos $middle >> $plot_base.ps\n")

    new_fp.write(b"if [ $montage != 0 ]\n")
    new_fp.write(b"then\n")
    new_fp.write(b"  title=${prefixes[$index]}\" \"$title\n")
    new_fp.write(b"  title_format=\"0 0 30 0 0 ML\" #Left-justify so montage titles are uniform.\n")
    new_fp.write(b"  title_x=$(bc <<< \"scale=5; $blurb_x-0.1\")\n")
    new_fp.write(b"else\n")
    new_fp.write(b"  title_x=$(bc <<< \"scale=5; $map_x+$map_width/2\")\n")
    new_fp.write(b"fi\n")
    new_fp.write(b"title_pos=\" -Xa${title_x}c -Ya${title_y}c \"\n")
    new_fp.write(b"echo $title_format $title | gmt pstext -N $title_pos $misc_range $middle >> $plot_base.ps\n")

    # Draw color scale with units printed above.
    write_gmt_colorscale(options, new_fp, 0)

    # Print blurb about data range, or masked amplitudes for phase plots.
    if len(results["rgb"]) == 3 and len(results["latlon"]["outputs"]) == 3:
        plot_options["blurb_disabled"] = 1
    if plot_options["blurb_disabled"]:
        new_fp.write(b"blurb_contents=\"\"\n")

    new_fp.write(b"echo $blurb_format $blurb_contents | gmt pstext -N $blurb_pos $misc_range $middle >> $plot_base.ps\n")

    blurb2_written = 0
    if len(results["latlon"]["outputs"]) == len(results["error_bars"]) and results["latlon"]["outputs"]:
        if results["error_bars"][i] > 0.0:
            blurb2_written = 1
            new_fp.write(f"blurb2_contents=\"Error bar: {results['error_bars'][i]:.1f} $scale_units\"\n".encode())

    if not blurb2_written:
        new_fp.write(b"blurb2_contents=\"\" #Error bar: N/A $scale_units\n")

    new_fp.write(b"echo $blurb_format $blurb2_contents | gmt pstext -N $blurb2_pos $misc_range $end >> $plot_base.ps\n")


def write_gmt_scripts(options: Options) -> None:
    """Generate GMT shell scripts for creating map plots.

    Args:
        options: Options object (plot_options, grid, results).
    """
    plot_options = options.plot_options
    grid         = options.grid
    results      = options.results

    new_file = Path(plot_options["outputfolder"]) / "create_plots.sh"
    try:
        new_fp = new_file.open("wb")
    except IOError:
        logging.error("The create_plots.sh GMT script couldn't be created.")

    new_fp.write(b"#!/bin/bash\n")
    new_fp.write(b"#set -x #Uncomment to echo these commands.\n")

    flip_file = Path(plot_options["outputfolder"]) / "flip_backgrounds.sh"
    try:
        flip_fp = flip_file.open("wb")
    except IOError:
        logging.error("The flip_backgrounds.sh script couldn't be created.")

    flip_fp.write(b"#!/bin/bash\n")
    flip_fp.write(b"set -x\n")

    trim_file = Path(plot_options["outputfolder"]) / "trim.sh"
    try:
        trim_fp = trim_file.open("wb")
    except IOError:
        logging.error("The trim.sh script couldn't be created.")

    trim_fp.write(b"#!/bin/bash\n")
    trim_fp.write(b"set -x\n")

    if 1:  #len(results['rgb']) == 3 and len(results['latlon']['outputs']) == 3:
        just_the_filenames = plot_options["just_the_filenames"][:1]

    for i in range(len(just_the_filenames)):
        if i == 0:
            write_gmt_defs(new_fp)
            new_fp.write(f"color_scheme={plot_options['color_scheme']} #1/2=white/black background\n".encode())
            new_fp.write(f"montage={plot_options['montage']} #1=left-justify titles, add (a),(b), run montage.sh.\n".encode())
            new_fp.write(b"prefixes=('(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' '(j)' '(k)' '(l)' '(m)' '(n)' '(o)' '(p)' '(q)' '(r)' '(s)' '(t)' '(u)' '(v)' '(w)' '(x)' '(y)' '(z)')\n")
            new_fp.write(b"index=-1 #Increments on each map, accesses prefixes above for montage.\n")
            new_fp.write(b'png_options=" -P -Tg " #PDF default: -E720, else 300 dpi.\n')
            new_fp.write(b"if [ $montage != 0 ]\nthen\n  png_options=\" -A\"$png_options\nfi\n")
            new_fp.write(b"#Force off-white(dark gray) fore(back)ground color because\n#flip_backgrounds.sh can change the maps' text from\n#black to white, and their backgrounds from white to black.\n")
            new_fp.write(b"gmt set COLOR_BACKGROUND=2/2/2 COLOR_FOREGROUND=253/253/253\n")
            new_fp.write(f"digits={plot_options['scale_digits']}\n".encode())
            new_fp.write(b"gmt set D_FORMAT=%.${digits}f\n")

        s = f"{plot_options['output_base']}_{i + 1:04d}"
        new_fp.write(b"#######################################################\n")
        if 1:  #len(results['rgb']) == 3 and len(results['latlon']['outputs']) == 3:
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
            flip_fp.write(f'{s}\n'.encode())
            trim_fp.write(f'{s}\n'.encode())
        else:
            flip_fp.write(f'{s}"\n'.encode())
            trim_fp.write(f'{s}"\n'.encode())

        write_gmt_map_data(options, new_fp, results["titles"][i], i)

        new_fp.write(b"convert $png_options $plot_base.ps #Convert PS to PNG format.\n")
        new_fp.write(b"#convert -P -Tf $plot_base.ps #Convert PS to PDF, if uncommented.\n")
        new_fp.write(b"rm -f $plot_base.ps\n")

        if len(results["rgb"]) != 3:
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
    extra_file = Path(plot_options["outputfolder"]) / "animate.sh"
    try:
        extra_fp = extra_file.open("wb")
    except IOError:
        logging.error("The animate.sh script couldn't be created.")

    extra_fp.write(b"#!/bin/bash\n")
    extra_fp.write(b"set -x\n")
    extra_fp.write(b"delay=100 #delay in hundredths of a second.\n")
    extra_fp.write(b"#size=\"640x480\"\n")
    extra_fp.write(b"#size=\"800x600\"\n")
    extra_fp.write(b"size=\"1024x768\"\n")
    extra_fp.write(f"output_base=\"{plot_options['output_base']}\"\n".encode())
    extra_fp.write(b"#Imagemagick can also output .mng (animated PNG, not well-supported), but ffmpeg is needed as a delegate for .mp4.\n")
    extra_fp.write(b"#convert -verbose -delay $delay -loop 0 $output_base* -resize $size animation.gif\n")
    extra_fp.write(b"#Or ffmpeg can output .mp4 directly.\n")
    extra_fp.write(b"ffmpeg -f image2 -i $output_base%d.png animation.mp4\n")
    extra_fp.close()

    # Create montage script
    extra_file = Path(plot_options["outputfolder"]) / "montage.sh"
    try:
        extra_fp = extra_file.open("wb")
    except IOError:
        logging.error("The montage.sh script couldn't be created.")

    extra_fp.write(b"#!/bin/bash\n")
    extra_fp.write(b"set -x\n")
    extra_fp.write(f"output_base=\"{plot_options['output_base']}\"\n".encode())
    extra_fp.write(b"montage $output_base* -geometry +2+2 montage.png\n")
    extra_fp.close()


def run_gmt_scripts(options: Options) -> None:
    """Set up and run GMT scripts via Docker.

    Args:
        options: Options object (output_folder).
    """
    docker_internal_folder   = "/home/jovyan"
    docker_internal_filename = "docker_internal.sh"
    docker_internal_script   = ("#!/bin/bash\n"
                                "gmt grdmix red.nc green.nc blue.nc -C -Gcombined.tif:GTiff\n"
                                "echo '----- Environment Variables -----'\n"
                                "echo 'PATH: '\n"
                                "echo $PATH\n"
                                "echo 'LD_LIBRARY_PATH: '\n"
                                "echo $LD_LIBRARY_PATH\n"
                                "echo '----- Conda Environments -----'\n"
                                "/opt/conda/bin/conda env list\n"
                                "echo '----- Working Directory -----'\n"
                                "pwd\n"
                                "cp create_plots.sh Zbackup_create_plots.sh\n"
                                "cp projections.sh Zbackup_projections.sh\n"
                                "./create_plots.sh\n")

    docker_command = (f"docker run -it -e TZ=America/Los_Angeles "
                      f"-v .:{docker_internal_folder} "
                      f"-w {docker_internal_folder} "
                      f"grace/testing-bpr-grace2 /bin/bash")

    docker_external_filename = "docker_external.bat"
    docker_external_script   = f"\n    {docker_command}\n    "

    # Write Docker internal script
    internal_path = options.output_folder / docker_internal_filename
    with internal_path.open("wb") as file:
        file.write(docker_internal_script.encode())

    # Write Docker external script
    external_path = options.output_folder / docker_external_filename
    with external_path.open("w") as file:
        file.write(docker_external_script)

    # Change permissions of the Docker internal script
    os.chmod(internal_path, 0o755)

    # Run Docker external script
    logging.info(f"Opening an interactive shell. Run {docker_external_filename}, then run {docker_internal_filename}")
    os.chdir(options.output_folder)
    os.system("cmd")


# ---------------------------------------------------------------------------
# Extraction helpers (called from main)
# ---------------------------------------------------------------------------

def zip_script(options: Options) -> None:
    """Zip the current script into the output folder.

    Args:
        options: Options object (output_folder).
    """
    current_script_path = Path(__file__).resolve()
    zip_file_path       = options.output_folder / (current_script_path.name + ".zip")
    with zipfile.ZipFile(zip_file_path, "w") as zipf:
        zipf.write(current_script_path, arcname=current_script_path.name)
    logging.info(f"Successfully zipped {os.fspath(current_script_path)} to {os.fspath(zip_file_path)}")


def init_plot_options(options: Options) -> None:
    """Initialize the plot_options dict on options.

    Args:
        options: Options object (output_folder, dpi).
    """
    options.plot_options = {
        "outputfolder"      : os.fspath(options.output_folder) + os.sep,
        "just_the_filenames" : ["r.nc", "g.nc", "b.nc"],
        "output_base"       : "map_parameter",
        "projection"        : 1,   # 1 = Robinson projection
        "plot_mascons"      : 0,
        "coastlines"        : 1,   # 1:coast, 2:coast+InSAR
        "blurb_disabled"    : 1,
        "montage"           : 0,   # 1=left-justify titles, add (a),(b), run montage.sh
        "region"            : "global",
        "frame"             : "a",
        "color_scheme"      : 2,   # 1/2=white/black background
        "land_color"        : "white",
        "sea_color"         : "blue",
        "scale_digits"      : 2,
        "show_fig"          : False,
        "save_fig"          : True,
        "figsize"           : (10, 10),
        "dpi"               : options.dpi,
        "linewidth"         : 2.0,
        "linestyle"         : "-",
        "color"             : "black",
        "marker"            : "o",
        "markersize"        : 5,
        "markerfacecolor"   : "blue",
        "markeredgewidth"   : 1.5,
        "markeredgecolor"   : "black",
        "font_size"         : 14,
        "font_weight"       : "bold",
        "x_label"           : "X-axis",
        "y_label"           : "Y-axis",
        "title"             : "My Plot",
        "grid"              : True,
        "grid_linestyle"    : "--",
        "grid_linewidth"    : 0.5,
        "grid_alpha"        : 0.7,
    }


def init_grid(options: Options) -> None:
    """Initialize the grid dict on options.

    Args:
        options: Options object.
    """
    options.grid = {
        "add_on" : {
            "write_gmt_defs"             : None,
            "just_the_filenames"         : [],
            "color_scheme"               : 1,   # 1/2=white/black background
            "montage"                    : 0,   # 1=left-justify titles, add (a),(b), run montage.sh
            "prefixes"                   : ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)",
                                            "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)",
                                            "(q)", "(r)", "(s)", "(t)", "(u)", "(v)", "(w)", "(x)",
                                            "(y)", "(z)"],
            "index"                      : -1,
            "png_options"                : " -P -Tg ",
            "digits"                     : 2,
            "output_base"                : "",
            "data_name"                  : "",
            "write_gmt_map_data"         : None,
            "convert"                    : None,
            "finish_flip_backgrounds"    : None,
            "finish_trim"                : None,
            "outputfolder"               : "",
            "animate.sh"                 : None,
            "montage.sh"                 : None,
            "write_gmt_colorscale"       : None,
            "write_rgb_colorscale"       : None,
        }
    }


def init_results(options: Options) -> None:
    """Initialize the results dict on options.

    Args:
        options: Options object (grid, plot_options, output_folder).
    """
    options.results = {
        "max_widths"         : 2,
        "options"            : {"output_choice": 5},
        "latlon"             : {
            "lat"        : [10.0, 20.0, 30.0, 40.0, 50.0],
            "lon"        : [10.0, 20.0, 30.0, 40.0, 50.0],
            "outputs"    : ["output1", "output2", "output3"],
            "mascon_lats" : [10.0, 20.0, 30.0],
        },
        "marker_lats"        : ["output1", "output2", "output3"],
        "marker_lons"        : ["output1", "output2", "output3"],
        "rgb"                : [1, 2, 3],
        "rgb_choice"         : 2,
        "units"              : ["unit1", "unit2", "unit3"],
        "maxlat"             : 90.0,
        "minlat"             : -90.0,
        "maxlon"             : 180.0,
        "minlon"             : -180.0,
        "error_bars"         : [0.1, 0.2, 0.3],
        "just_the_filenames" : ["example_filename"],
        "color_scheme"       : 1,
        "montage"            : 1,
        "scale_digits"       : 2,
        "output_base"        : "output",
        "titles"             : ["example_title"],
        "grid"               : options.grid,
        "plot_options"       : options.plot_options,
        "outputfolder"       : os.fspath(options.output_folder),
        "date"               : dt.datetime.now(),
        "config"             : {"x": 0.0, "y": 0.0, "width": 0.0, "height": 0.0},
        "xy"                 : {"x_values": [[0.0]], "y_values": [[0.0]]},
        "xyz"                : {"x_values": [], "y_values": [], "z_values": []},
        "min_max"            : {"x": [0.0, 0.0], "y": [0.0, 0.0], "z": [0.0, 0.0]},
        "misc"               : {
            "scale_format" : "",
            "overflow"     : "",
            "scale_pos"    : "",
            "units_format" : "",
            "units_pos"    : "",
            "misc_range"   : "",
            "start"        : "",
            "middle"       : "",
            "end"          : "",
            "scale_width"  : 0.0,
            "scale_x"      : 0.0,
            "scale_y"      : 0.0,
            "scale_length" : 0.0,
            "plot_base"    : "",
        },
    }


def configure_keys_for_input(options: Options) -> None:
    """Set data key names and per-dataset defaults based on input_choice.

    Args:
        options: Options object. Mutates x_key, y_key, lat_key, lon_key,
                 min_period, max_period, xskip.
    """
    if options.input_choice == "SSHA":
        options.x_key      = "Time"
        options.y_key      = "SLA"
        options.lat_key    = "Latitude"
        options.lon_key    = "Longitude"
        options.min_period = 30.0
        options.max_period = 60.0
        options.xskip      = 6
    elif options.input_choice == "Argo":
        options.x_key      = "time"
        options.y_key      = "ohc_2d_lt_700m"
        options.lat_key    = "latitude"
        options.lon_key    = "longitude"
        options.min_period = 55 * 7.0
        options.max_period = 150 * 7.0
        options.xskip      = 1
    elif options.input_choice == "MUR_SST":
        options.x_key      = "time"
        options.y_key      = "sst_anomaly"
        options.lat_key    = "lat"
        options.lon_key    = "lon"
        options.min_period = 2 * 7.0
        options.max_period = 24 * 7.0
        options.xskip      = 40
    elif options.input_choice == "AQUA_MODIS":
        options.x_key      = "time"
        options.y_key      = "sst_anomaly"
        options.lat_key    = "lat"
        options.lon_key    = "lon"
        options.min_period = 2 * 7.0
        options.max_period = 24 * 7.0
        options.xskip      = 40
    elif options.input_choice == "GRACE":
        options.x_key      = "time"
        options.y_key      = "lwe_thickness"
        options.lat_key    = "lat"
        options.lon_key    = "lon"
        options.min_period = 24 * 7.0
        options.max_period = 70 * 7.0
        options.xskip      = 1
    else:
        logging.error(f"DID NOT RECOGNIZE {options.input_choice = }")


def load_input_data(options: Options) -> xr.Dataset:
    """Load input data based on the configured input_choice.

    Args:
        options: Options object (input_choice, and all folder paths).

    Returns:
        Loaded xarray Dataset.
    """
    if options.input_choice == "SSHA":
        input_data = load_ssha_files(options, tskip=1)
    elif options.input_choice == "Argo":
        input_data = load_argo_file(options)
        input_data["ohc_2d_lt_700m"] = input_data["ohc_2d"] - input_data["ohc_2d_700m"]
    elif options.input_choice == "MUR_SST":
        input_data = load_mur_sst_files(options, tskip=1)
    elif options.input_choice == "AQUA_MODIS":
        input_data = load_aqua_modis_files(options, tskip=1)
        logging.info(input_data)
    elif options.input_choice == "GRACE":
        input_data = load_grace_file(options)
    else:
        logging.error(f"DID NOT RECOGNIZE {options.input_choice = }")
        raise ValueError(f"Unknown input_choice: {options.input_choice}")
    return input_data


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
