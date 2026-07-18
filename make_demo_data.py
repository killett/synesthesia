"""Generate synthetic weekly sea-surface grids for the quickstart.

Two regions oscillate at different periods (60 and 120 days), so the
output map shows them as two distinct hues.
Run from the repo root: .venv/bin/python make_demo_data.py
"""
from pathlib import Path

import numpy as np
import xarray as xr

rng = np.random.default_rng(42)
lat = np.arange(-88.0, 90.0, 4.0)
lon = np.arange(0.0, 360.0, 4.0)
weeks = np.arange(
    np.datetime64("2024-01-03"),
    np.datetime64("2026-01-01"),
    np.timedelta64(7, "D"),
)

la, lo = np.meshgrid(lat, lon, indexing="ij")
box_60d = (la > 20) & (la < 60) & (lo > 120) & (lo < 200)
box_120d = (la < -10) & (la > -50) & (lo > 240) & (lo < 320)

for i, t in enumerate(weeks):
    days = float((t - weeks[0]) / np.timedelta64(1, "D"))
    field = 0.02 * rng.standard_normal(la.shape)
    field += np.where(box_60d, 0.5 * np.sin(2 * np.pi * days / 60.0), 0.0)
    field += np.where(box_120d, 0.5 * np.sin(2 * np.pi * days / 120.0), 0.0)
    ds = xr.Dataset(
        {
            "ssha": (("latitude", "longitude"), field.astype(np.float32)),
            "time": ((), t),  # scalar variable, not a coord — the loader promotes it
        },
        coords={"latitude": lat, "longitude": lon},
    )
    out = (
        Path("input_timeseries/simple_grids")
        / str(t.astype("datetime64[Y]"))
        / f"week_{i:03d}.nc"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out)

print(f"wrote {len(weeks)} files")
