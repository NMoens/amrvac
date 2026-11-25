#!/usr/bin/env python3
"""
Compare the 1D resolution-scan runs (resscan/nXXXX_NNNN.dat).

  python3 resscan/plot_scan.py [--snap 100] [--res 0320 0640 1280] [-o scan_profiles.png]

Plots density, radial velocity, temperatures (T_gas solid, T_rad dashed) and
radiative luminosity against r for one snapshot of each resolution: the full
domain and a zoom on the photosphere / wind base (cells marked), and prints
a table of photospheric radius (tau_F = 2/3), mass-loss rate and velocities.
Reads the .dat files directly (double precision, conservative variables).
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'm1_Ostar_2D'))
from dat_lift import read_snapshot, map_to_level  # noqa: E402

# code units of the Star8 setup (printed at the start of every run)
UNIT_L, UNIT_RHO, UNIT_V, UNIT_P = 7.6527e11, 3.8e-6, 8.2343181e6, 2.5765518e8
MSUN, YEAR, A_RAD = 1.98847e33, 3.15576e7, 7.5657e-15
MDOT = 4 * np.pi * UNIT_RHO * UNIT_V * UNIT_L ** 2 / MSUN * YEAR   # code -> Msun/yr
# ordinal one-hue ramp (light = coarse, dark = fine)
COLORS = {'0320': '#86b6ef', '0640': '#3987e5', '1280': '#1c5cab', '2560': '#0d366b'}


def profile(fn):
    h, _, b = read_snapshot(fn)
    U = map_to_level(h, b, 1, False)
    n = h['w_names']
    nx = int(h['domain_nx'][0])
    e = h['xmin'][0] * (h['xmax'][0] / h['xmin'][0]) ** (np.arange(nx + 1) / nx)
    r = 0.5 * (e[1:] + e[:-1])
    rho = U[:, n.index('rho')]
    v = U[:, n.index('m1')] / rho
    out = dict(t=h['time'], r=r, dr=np.diff(e), rho=rho, v=v,
               Trad=(U[:, n.index('Er')] * UNIT_P / A_RAD) ** 0.25,
               mdot=MDOT * r ** 2 * rho * v)
    # diagnostics stored by the run (computed at output time)
    for k in ('Tgas', 'lum'):
        if k in n:
            out[k] = U[:, n.index(k)]
    if 'kappa_F' in n:
        dtau = rho * U[:, n.index('kappa_F')] * out['dr']
        out['tau'] = np.cumsum(dtau[::-1])[::-1]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--snap', type=int, default=100, help='snapshot index (100 = t=1)')
    ap.add_argument('--res', nargs='+', default=['0320', '0640', '1280'])
    ap.add_argument('--zoom', nargs=2, type=float, default=[1.15, 1.45], help='r range of the zoom')
    ap.add_argument('-o', '--output', default=os.path.join(HERE, 'scan_profiles.png'))
    args = ap.parse_args()

    P = {}
    for res in args.res:
        fn = os.path.join(HERE, 'n%s_%04d.dat' % (res, args.snap))
        if not os.path.exists(fn):
            print('missing %s, skipped' % fn)
            continue
        P[res] = profile(fn)

    # rows: (y data function, y label, log-scale y). Temperatures: T_gas solid,
    # T_rad dashed; colour always encodes the resolution
    rows = [
        (lambda p: [('', np.log10(p['rho']))], 'log$_{10}$ $\\rho$ (code units)', False),
        (lambda p: [('', p['v'] * UNIT_V / 1e5)], '$v_r$ (km/s)', False),
        (lambda p: [('T$_{gas}$', p['Tgas'] / 1e3), ('T$_{rad}$', p['Trad'] / 1e3)],
         'T (kK)', True),
        (lambda p: [('', p['lum'])], 'L$_{rad}$ / L$_*$', False),
    ]
    t0 = next(iter(P.values()))['t']
    fig, axes = plt.subplots(len(rows), 2, figsize=(12, 3.4 * len(rows)), sharex='col')
    for irow, (getter, ylab, logy) in enumerate(rows):
        ax1, ax2 = axes[irow]
        for res, p in P.items():
            col = COLORS.get(res, 'k')
            m = (p['r'] > args.zoom[0]) & (p['r'] < args.zoom[1])
            for k, (name, y) in enumerate(getter(p)):
                ls = '-' if k == 0 else '--'
                lab = '%d cells' % int(res) if k == 0 else None
                ax1.plot(p['r'], y, color=col, lw=2 if k == 0 else 1.5, ls=ls, label=lab)
                ax2.plot(p['r'][m], y[m], color=col, lw=2 if k == 0 else 1.5, ls=ls,
                         marker='o' if k == 0 else None, ms=4, mec='white', mew=0.8, label=lab)
        for ax in (ax1, ax2):
            ax.set_ylabel(ylab)
            if logy:
                ax.set_yscale('log')
            ax.grid(alpha=0.25, lw=0.6)
            for s in ('top', 'right'):
                ax.spines[s].set_visible(False)
        ax2.set_xlim(*args.zoom)
        # zoom y-range: data inside the zoom only
        ys = np.concatenate([y[(p['r'] > args.zoom[0]) & (p['r'] < args.zoom[1])]
                             for p in P.values() for _, y in getter(p)])
        lo, hi = np.nanmin(ys), np.nanmax(ys)
        if logy:
            ax2.set_ylim(lo / 1.15, hi * 1.15)
        else:
            pad = 0.05 * (hi - lo if hi > lo else 1.0)
            ax2.set_ylim(lo - pad, hi + pad)
    axes[0, 0].set_title('Whole domain, t = %.2f' % t0, fontsize=11)
    axes[0, 1].set_title('Photosphere and wind base (one marker per cell)', fontsize=11)
    axes[0, 0].legend(frameon=False)
    # explain the line styles in the temperature row
    from matplotlib.lines import Line2D
    axes[2, 0].legend(handles=[Line2D([], [], color='0.35', lw=2, ls='-', label='T$_{gas}$'),
                               Line2D([], [], color='0.35', lw=1.5, ls='--', label='T$_{rad}$')],
                      frameon=False)
    for ax in axes[-1]:
        ax.set_xlabel('r / R$_*$')
    fig.tight_layout()
    fig.savefig(args.output, dpi=110)
    print('wrote', args.output)

    # table
    print('\n%-6s %6s %9s %9s %9s %9s %10s %10s' % ('cells', 't', 'R_ph', 'dr(R_ph)',
          'Mdot1.5', 'Mdot4.5', 'v(1.5)', 'v(4.5)'))
    for res, p in P.items():
        r = p['r']
        rph = np.interp(2 / 3, p['tau'][::-1], r[::-1]) if 'tau' in p else np.nan
        k = np.argmin(abs(r - rph))
        at = lambda q, x: np.interp(x, r, p[q])
        print('%-6s %6.3f %9.4f %9.5f %9.2e %9.2e %10.1f %10.1f' % (
            res, p['t'], rph, p['dr'][k], at('mdot', 1.5), at('mdot', 4.5),
            at('v', 1.5) * UNIT_V / 1e5, at('v', 4.5) * UNIT_V / 1e5))
    print('(R_ph at tau_F = 2/3; Mdot in Msun/yr; v in km/s)')


if __name__ == '__main__':
    main()
