# synesthesia

Synesthesia renders the frequency content of gridded earth-science time
series as colour. For each pixel of a stack of NetCDF snapshots it computes an
NFFT power spectrum over a chosen band of oscillation periods, then projects
that spectrum through the CIE 1931 colour-matching functions as if the
temporal frequencies were visible light — so different oscillation periods
land at different hues. The result is an sRGB map, rendered with
matplotlib/cartopy, in which a region that wobbles with a 60-day period looks a
different colour from one that wobbles with a 120-day period. The core is
`timeseries2color.py`, a Python port of the Hughes & Williams (2010) approach.

## Table of contents

- [How it works](#how-it-works)
- [Quickstart](#quickstart)
- [Installation](#installation)
- [Command cheatsheet](#command-cheatsheet)
- [Repository tour](#repository-tour)
- [Troubleshooting](#troubleshooting)
- [Provenance & license](#provenance--license)

## How it works

1. **Load a gridded time series.** One NetCDF snapshot per timestamp is read
   and stacked along time. Supported sources are selected with
   `--input-choice`: `SSHA` (sea-surface height anomaly), `GRACE`, `MUR_SST`,
   `AQUA_MODIS` (Aqua MODIS chlorophyll), `Argo` (ocean heat content), and
   `simple_grids` (a generic recursive folder of NetCDF snapshots, used by the
   quickstart below).
2. **Compute a per-pixel spectrum.** For every grid cell, an NFFT
   (non-uniform FFT) power spectrum is taken over the requested band of
   oscillation periods (`--min-period` to `--max-period`, in days). The number
   of timestamps is reduced to an even count first, because NFFT requires it —
   you will see a log line reporting this.
3. **Map the spectrum to colour.** Each pixel's power spectrum is projected
   onto the CIE 1931 colour-matching functions (the CSV tables in
   `CIE_1931/`), exactly as if the temporal spectrum were a spectrum of
   light. This yields CIE XYZ, which is converted to sRGB.
4. **Render the map.** The sRGB field is drawn as a world map with
   matplotlib and cartopy (with a configurable projection, colours, coastlines
   and so on), and — unless `--no-gmt` is passed — an accompanying set of GMT
   plotting scripts is generated alongside the PNG.

## Quickstart

This walks through generating synthetic demo data and rendering it, end to
end, with no external datasets. It assumes you have completed
[Installation](#installation) first (the venv at `.venv`, `nfft` installed).
All commands below are run from the repo root.

**1. Generate the demo data.** The repo ships `make_demo_data.py`, which
writes 105 weekly 45×90 sea-surface grids in which two rectangular regions
oscillate at 60-day and 120-day periods:

```bash
.venv/bin/python make_demo_data.py
```

It prints `wrote 105 files`, written under `input_timeseries/simple_grids/`.

**2. Render the map:**

```bash
.venv/bin/python timeseries2color.py "quickstart demo" \
    --input-choice simple_grids --xskip 1 --no-gmt \
    --min-period 30 --max-period 160
```

`--xskip 1` is load-bearing: the default decimates the grid heavily (see
[Troubleshooting](#troubleshooting)). `--no-gmt` skips the GMT script pipeline,
which needs system GMT binaries. `--min-period 30 --max-period 160` brackets
the demo's 60-day and 120-day signals.

**What you should see.** The run takes about half a minute on a typical
machine. Output lands in a
timestamped folder, `output/<YYYYMMDD-HHMMSS> - quickstart demo/`, containing
`simple_grids_map_matplotlib.png`. On a dark background the map shows two
solid-coloured rectangles: a northern box (the 60-day region) rendered blue,
and a southern box (the 120-day region) rendered red/orange. The two boxes are
distinct hues because their oscillation periods differ — that is the whole
point of the tool.

## Installation

Instructions verified on linux-64 with Python 3.12.13, using
[uv](https://docs.astral.sh/uv/) to manage the environment; the kerchunk and
GMT paths below are optional and were not verified.

**1. Create the virtual environment.** If `uv` cannot find a local Python
3.12 it downloads one automatically:

```bash
uv venv --python 3.12 .venv
```

**2. Install the pinned dependencies:**

```bash
uv pip install --python .venv/bin/python \
    "numpy==2.4.4" "scipy==1.17.1" "pandas==3.0.1" \
    "xarray==2026.2.0" "h5py==3.16.0" "netCDF4==1.7.4" \
    "dask==2026.3.0" "statsmodels==0.14.6" "matplotlib==3.10.8" \
    "cartopy>=0.24" "colour-science==0.4.7" "tqdm>=4.66"
```

**3. Install nfft** (from PyPI; builds from source, which is expected):

```bash
uv pip install --python .venv/bin/python nfft
```

**4. Verify** the tool imports and reports its version:

```bash
.venv/bin/python timeseries2color.py --version
```

This prints `timeseries2color.py 0.2.1`.

**Optional — kerchunk.** The `--use-kerchunk` loader path needs the `kerchunk`
package, which is not in the pinned list above. Install it only if you intend
to use that flag:

```bash
# optional, only for --use-kerchunk:
# uv pip install --python .venv/bin/python kerchunk
```

**Optional — GMT system packages.** GMT script generation (everything except
`--no-gmt`) needs GMT and a few companion tools installed at the system level,
only required when you do **not** pass `--no-gmt`:

```bash
# optional, only needed without --no-gmt:
# sudo apt-get install gmt gmt-dcw gmt-gshhg ghostscript bc imagemagick
```

## Command cheatsheet

Every flag below comes straight from `timeseries2color.py --help` (version
0.2.1). Defaults shown as "per-dataset" are chosen automatically based on
`--input-choice`.

### Input selection

| Flag | Values / default | What it does |
|------|------------------|--------------|
| `description` (positional) | free text | Free-text run label (quote if it contains spaces); appears in the output folder name. |
| `--input-choice` | `{SSHA, simple_grids, Argo, MUR_SST, AQUA_MODIS, GRACE}` (default: `AQUA_MODIS`) | Data source. |
| `--lat-var` | name | Override auto-detected latitude variable name. |
| `--lon-var` | name | Override auto-detected longitude variable name. |
| `--time-var` | name | Override auto-detected time variable name. |
| `--use-kerchunk` | flag | Load multi-file datasets through a kerchunk virtual-zarr reference. Faster for large xskip values; builds a sidecar `.kerchunk.json` on first run. |
| `--kerchunk-path` | path (default: `<data folder>/.kerchunk.json`) | Override the kerchunk reference file path. |

### Analysis

| Flag | Values / default | What it does |
|------|------------------|--------------|
| `--xskip` | int (default: per-dataset) | Skip every N points in lat/lon. |
| `--min-period` | days (default: per-dataset) | Minimum period in days. |
| `--max-period` | days (default: per-dataset) | Maximum period in days. |
| `--chunk-bytes` | bytes (default: 200 MB) | Target bytes per dask time chunk. |
| `--block-bytes` | bytes (default: 500 MB) | Target bytes per loop lat/lon block. |
| `--max-graph-chunks` | int (default: 1,000,000) | Hard cap on total dask chunks across all opened files. Lower if `open_mfdataset` OOMs; raise for finer-grained spatial chunks. |
| `--lat-chunk` | int | Override auto-derived lat chunk size for `open_mfdataset`. Must be paired with `--lon-chunk`. |
| `--lon-chunk` | int | Override auto-derived lon chunk size for `open_mfdataset`. Must be paired with `--lat-chunk`. |
| `--fit-terms` | one or more of `{constant, trend, accel, annual, semiannual}` (default: all five) | Detrending terms. |

### Output & plotting

| Flag | Values / default | What it does |
|------|------------------|--------------|
| `--dpi` | int (default: 300) | DPI for output images. |
| `-d`, `--debug` | flag | Enable debug logging. |
| `--no-gmt` | flag | Skip GMT plot generation. |
| `--projection` | `{1, 2, 3, 4}` | Map projection: 1=Robinson, 2=Winkel Tripel, 3=Mollweide, 4=Miller. |
| `--use-old-funcs` | flag | Use Hughes & Williams 2010 functions instead of colour-science. (This is currently the default path; there is no flag to enable the colour-science path, so this flag is effectively a no-op.) |
| `--land-color` | colour (default: #333333) | Land color for matplotlib maps. |
| `--ocean-color` | colour (default: black) | Ocean color for matplotlib maps. |
| `--light-mode` | flag | Use light/white background instead of dark. |
| `--borders` | flag | Show political boundaries on maps. |
| `--rivers` | flag | Show major rivers on maps. |
| `--coastline-resolution` | `{110m, 50m}` (default: 110m) | Coastline resolution. |

`-h`/`--help` and `-v`/`--version` print help and the version and exit.

## Repository tour

- `timeseries2color.py` — the tool (~3300 lines); everything the quickstart
  and cheatsheet describe.
- `make_demo_data.py` — generates the quickstart's synthetic demo data.
- `notation.sh`, `overflow.sh`, `projections.sh` — GMT plotting helpers
  consumed by the generated GMT scripts.
- `CIE_1931/` — the CIE 1931 colour-matching tables.
- `input_timeseries/` — root for input data; created on first use, nothing
  in it is tracked. The quickstart's generated demo data lands in
  `input_timeseries/simple_grids/`.

## Troubleshooting

**My map is almost blank / only a point or two of colour.** *Cause:*
`--xskip` defaults to a per-dataset value (100 for the generic path), which
decimates a small grid down to almost nothing. *Fix:* pass a smaller value; for
the quickstart's 45×90 grid use `--xskip 1`.

**A GMT command is not found, or the run fails during GMT script execution
when you did not pass `--no-gmt`.** *Cause:* by default the tool
generates and runs GMT plotting scripts, which require system GMT binaries that
are not present in every environment. *Fix:* either install the optional GMT
system packages listed under [Installation](#installation), or pass `--no-gmt`
to skip GMT entirely and produce only the matplotlib/cartopy PNG.

**`ModuleNotFoundError: No module named 'kerchunk'` when you pass
`--use-kerchunk`.** *Cause:* `kerchunk` is imported lazily only when that flag
is used, and it is not part of the pinned dependency list. *Fix:* install it
with `uv pip install --python .venv/bin/python kerchunk`.

*Not an error:* during `uv pip install` you may see a warning that uv "Failed
to hardlink files; falling back to full copy" — this is a performance note
(the cache and target directories are on different filesystems), not a failure,
and can be silenced with `export UV_LINK_MODE=copy` if desired.

## Provenance & license

The original C++ implementation of Hughes & Williams (2010) lives in the
[GAIA](https://github.com/killett/GAIA) repository.

Apache-2.0 — see [LICENSE](LICENSE).
