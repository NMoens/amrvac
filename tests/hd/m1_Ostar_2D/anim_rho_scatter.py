#!/usr/bin/env python3
"""
Animate a cell-by-cell scatter plot of log10(rho) (or another variable) versus
radius for a series of 2D spherical (r, theta) MPI-AMRVAC .vtu snapshots.

Every frame is one snapshot (one time); every point is one grid cell. The .vtu
files store Cartesian cell coordinates (x = r sin(theta), y = r cos(theta)),
so r and theta are reconstructed from the cell centres. Axis limits are fixed
over all frames so that the evolution can be followed.

Examples:
  python3 anim_rho_scatter.py output/Star8_*.vtu -o rho_r.gif
  python3 anim_rho_scatter.py output/Star8_*.vtu -o rho_r.gif --color none --every 4
  python3 anim_rho_scatter.py output/Star8_*.vtu -o Trad_r.gif --var Trad
  python3 anim_rho_scatter.py output/Star8_*.vtu -o rho_zoom.gif --rmax 2 --ymin -10
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def read_vtu(fname, var, need_level=False):
    """Return r, theta, value, theta cell width (or None) and time of one .vtu file."""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fname)
    reader.Update()
    grid = reader.GetOutput()

    cdata = grid.GetCellData()
    arr = cdata.GetArray(var)
    if arr is None:
        names = [cdata.GetArrayName(i) for i in range(cdata.GetNumberOfArrays())]
        sys.exit("Variable '%s' not found in %s. Available: %s" % (var, fname, ", ".join(names)))
    val = vtk_to_numpy(arr).astype(np.float64)

    centres = vtk.vtkCellCenters()
    centres.SetInputData(grid)
    centres.Update()
    xy = vtk_to_numpy(centres.GetOutput().GetPoints().GetData())
    r = np.hypot(xy[:, 0], xy[:, 1])
    theta = np.arctan2(xy[:, 0], xy[:, 1])

    dtheta = None
    if need_level:
        pts = vtk_to_numpy(grid.GetPoints().GetData()).astype(np.float64)
        conn = vtk_to_numpy(grid.GetCells().GetConnectivityArray()).reshape(-1, 4)
        tp = np.arctan2(pts[conn, 0], pts[conn, 1])
        dtheta = tp.max(axis=1) - tp.min(axis=1)

    tarr = grid.GetFieldData().GetArray("TIME")
    time = float(vtk_to_numpy(tarr)[0]) if tarr is not None else np.nan
    return r, theta, val, dtheta, time


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+", help=".vtu snapshots (sorted by name)")
    p.add_argument("-o", "--output", default="rho_r.gif", help="output gif (default rho_r.gif)")
    p.add_argument("--var", default="rho", help="cell variable to plot (default rho)")
    p.add_argument("--linear", action="store_true", help="plot the value itself instead of log10")
    p.add_argument("--color", choices=["theta", "level", "none"], default="theta",
                   help="colour points by theta/pi, by AMR level, or not at all (default theta)")
    p.add_argument("--every", type=int, default=1, help="plot every N-th cell (default 1)")
    p.add_argument("--every-frame", type=int, default=1, help="use every K-th snapshot (default 1)")
    p.add_argument("--rmin", type=float, help="lower radius limit")
    p.add_argument("--rmax", type=float, help="upper radius limit")
    p.add_argument("--ymin", type=float, help="lower y limit (in plotted units, i.e. log10 by default)")
    p.add_argument("--ymax", type=float, help="upper y limit")
    p.add_argument("--fps", type=float, default=4, help="frames per second (default 4)")
    p.add_argument("--dpi", type=int, default=100, help="resolution (default 100)")
    p.add_argument("--size", type=float, default=1.0, help="marker size (default 1)")
    args = p.parse_args()

    files = sorted(set(args.files))[::max(1, args.every_frame)]
    if not files:
        sys.exit("No input files")

    # Pass 1: read all frames once, cache the (subsampled) plot data
    frames = []
    for i, f in enumerate(files):
        r, th, val, dth, t = read_vtu(f, args.var, need_level=(args.color == "level"))
        sl = slice(None, None, max(1, args.every))
        r, th, val = r[sl], th[sl], val[sl]
        if args.linear:
            y = val
        else:
            with np.errstate(divide="ignore", invalid="ignore"):
                y = np.where(val > 0, np.log10(val), np.nan)
        if args.color == "theta":
            c = th / np.pi
        elif args.color == "level":
            dth = dth[sl]
            c = np.round(np.log2(dth.max() / dth)) + 1
        else:
            c = None
        ok = np.flatnonzero(np.isfinite(y))
        # Random drawing order, so no theta range is systematically hidden under another
        ok = ok[np.random.default_rng(i).permutation(ok.size)]
        frames.append(dict(r=r[ok].astype(np.float32), y=y[ok].astype(np.float32),
                           c=None if c is None else c[ok].astype(np.float32),
                           t=t, name=os.path.basename(f)))
        print("read %3d/%d  %s  t=%.4g  (%d points)" % (i + 1, len(files), f, t, ok.sum()))

    # Fixed axis limits over all frames
    allr = np.concatenate([fr["r"] for fr in frames])
    ally = np.concatenate([fr["y"] for fr in frames])
    rmin = args.rmin if args.rmin is not None else allr.min()
    rmax = args.rmax if args.rmax is not None else allr.max()
    ylo, yhi = np.percentile(ally, [0.1, 100.0])
    pad = 0.03 * (yhi - ylo)
    ymin = args.ymin if args.ymin is not None else ylo - pad
    ymax = args.ymax if args.ymax is not None else yhi + pad

    # Pass 2: animation
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_xlim(rmin, rmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("r")
    ax.set_ylabel(args.var if args.linear else "log10 " + args.var)
    ax.grid(alpha=0.3)

    f0 = frames[0]
    if args.color == "none":
        sc = ax.scatter(f0["r"], f0["y"], s=args.size, c="C0", lw=0, rasterized=True)
    else:
        allc = np.concatenate([fr["c"] for fr in frames])
        cmap = "turbo" if args.color == "theta" else "viridis"
        sc = ax.scatter(f0["r"], f0["y"], s=args.size, c=f0["c"], cmap=cmap,
                        vmin=allc.min(), vmax=allc.max(), lw=0, rasterized=True)
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label("theta / pi" if args.color == "theta" else "AMR level")
    # Lay out with a representative title so it is not clipped
    title = ax.set_title("%s   t = %.4g" % (frames[-1]["name"], frames[-1]["t"]))
    fig.tight_layout()

    def update(i):
        fr = frames[i]
        sc.set_offsets(np.column_stack([fr["r"], fr["y"]]))
        if fr["c"] is not None:
            sc.set_array(fr["c"])
        title.set_text("%s   t = %.4g" % (fr["name"], fr["t"]))
        return sc, title

    # Render every frame to RGB, then write the gif with one shared palette and
    # no dithering (matplotlib's PillowWriter quantises each frame separately,
    # which gives washed-out, dithered frames)
    fig.set_dpi(args.dpi)
    images = []
    for i in range(len(frames)):
        update(i)
        fig.canvas.draw()
        images.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()))
    plt.close(fig)

    palette_src = Image.new("RGB", (images[0].width, images[0].height * len(images)))
    for k, im in enumerate(images):
        palette_src.paste(im, (0, k * im.height))
    pal = palette_src.quantize(colors=256, method=Image.Quantize.MEDIANCUT)
    gif = [im.quantize(palette=pal, dither=Image.Dither.NONE) for im in images]
    gif[0].save(args.output, save_all=True, append_images=gif[1:],
                duration=int(round(1000.0 / args.fps)), loop=0, optimize=False)
    print("wrote %s (%d frames, %.1f MB)" % (args.output, len(frames), os.path.getsize(args.output) / 1e6))


if __name__ == "__main__":
    main()
