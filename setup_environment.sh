#!/usr/bin/env bash
#
# setup_environment.sh — Install uv, create a Python 3.12 venv, and install
# all dependencies for pyGAIA-copy on Ubuntu 24.04 LTS.
#
# Run from the pyGAIA-copy/ directory:
#     chmod +x setup_environment.sh
#     ./setup_environment.sh
#
# This script uses uv instead of pip, which sidesteps the PEP 668
# "externally managed environment" restriction on Ubuntu 24.04.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# ---------- 1. Install uv (standalone installer, no pip needed) ----------

if command -v uv &>/dev/null; then
    echo "uv is already installed: $(uv --version)"
else
    echo "Installing uv..."
    curl -LsSf https://astral.sh/uv/install.sh | sh

    # The installer puts uv in ~/.local/bin (or ~/.cargo/bin).
    # Make it available in this session.
    export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"

    if ! command -v uv &>/dev/null; then
        echo "ERROR: uv not found on PATH after install. Check the installer output above."
        exit 1
    fi
    echo "Installed uv: $(uv --version)"
fi

# ---------- 2. Create virtual environment with Python 3.12 ----------

if [ -d ".venv" ]; then
    echo ".venv already exists — skipping creation."
    echo "  (Delete it first if you want a fresh environment: rm -rf .venv)"
else
    echo "Creating .venv with Python 3.12..."
    uv venv --python 3.12 .venv
fi

# ---------- 3. Install pinned PyPI packages ----------

echo "Installing Python packages..."
uv pip install --python .venv/bin/python \
    "numpy==2.4.4" \
    "scipy==1.17.1" \
    "pandas==3.0.1" \
    "xarray==2026.2.0" \
    "h5py==3.16.0" \
    "netCDF4==1.7.4" \
    "dask==2026.3.0" \
    "statsmodels==0.14.6" \
    "matplotlib==3.10.8" \
    "cartopy>=0.24" \
    "colour-science==0.4.7" \
    "tqdm>=4.66"

# ---------- 4. Install nfft from local source ----------

if [ -d "nfft-source/nfft-master" ]; then
    echo "Installing nfft from local source..."
    uv pip install --python .venv/bin/python ./nfft-source/nfft-master/
else
    echo "WARNING: nfft-source/nfft-master/ not found — skipping nfft install."
    echo "  You can install it later with:"
    echo "    uv pip install --python .venv/bin/python ./nfft-source/nfft-master/"
fi

# ---------- 5. Install system dependencies for GMT pipeline ----------

echo ""
echo "Installing system packages (requires sudo)..."
sudo apt-get update -qq
sudo apt-get install -y gmt gmt-dcw gmt-gshhg ghostscript bc imagemagick

# ---------- Done ----------

echo ""
echo "=========================================="
echo "  Environment setup complete."
echo ""
echo "  Activate with:"
echo "    source .venv/bin/activate"
echo ""
echo "  Test run:"
echo "    python functions49.py --xskip 192"
echo "=========================================="
