#!/usr/bin/env python3
"""
Lift a spherical MPI-AMRVAC snapshot (.dat) to more dimensions, to restart a
higher-dimensional run from it:  1D -> 2D,  1D -> 3D,  2D -> 3D.

The solution is copied along every new dimension (theta for 1D sources, phi
for 3D targets); the new vector components (m2/F2 for 1D sources, m3/F3 for
2D sources without them) are zero. All output blocks are at level 1.

Kept dimensions (r, and theta for 2D sources)
  Taken from the source file: xprobmin/xprobmax and domain_nx * 2**(L-1),
  with L the level the source is mapped to (default: the finest level in the
  file; --level to choose). Coarser source blocks are copied, finer ones are
  volume-averaged with spherical weights (r^2 dr, sin(theta) dtheta). The
  radial grid is assumed to be amrvac's 'uni' stretched grid with the default
  qstretch (or uniform with --uniform-r). Block sizes may differ from the
  source (--bnx1, --bnx2 or the par files).
  NOTE: the target run must use the same radial stretching as the source run
  (stretch_dim(1), qstretch_baselevel); this is not stored in .dat files.

New dimensions (theta for 1D sources, phi for 3D targets)
  xprobmin/xprobmax (par-file units: multiples of 2*pi), domain_nx and
  block_nx from the target par file(s) given with --par (later files override
  earlier ones), and/or the command line (--theta-min/--theta-max/--nx2/--bnx2,
  --phi-min/--phi-max/--nx3/--bnx3), which override the par files.
  Periodicity of a new dimension (stored in the .dat and checked at restart):
  --periodic-phi/--no-periodic-phi, else detected from typeboundary_min3/max3
  ('periodic') in the par files, else non-periodic.

Variables
  Default: rho, m1..m<ndir>, e, Er, F1..F<ndir> with ndir = target ndim
  (--ndir 3 for a 2.5D target), i.e. the M1 hd flux variables in amrvac order;
  they must be the first variables of the target run. Copied by name from the
  source; m2, m3, F2, F3 are zero if absent. --vars overrides the list.
  Extra variables (opacities, diagnostics) are recomputed by amrvac.

Outputs
  * the .dat file (format version 5), default <source base><N>d<NNNN>.dat
  * a .par file with restart_from_file and the full mesh, to be given LAST:
        ./amrvac -i <base>.par <star>.par <generated>.par

Examples
  # 1D -> 2D (theta grid from the 2D par files), start at t=0
  python3 dat_lift.py ../m1_Ostar_1D/resscan/n1280_0100.dat --ndim 2 --reset \\
          --par parfiles_hpc/parfiles/base.par parfiles_hpc/parfiles/Star8.par
  # 2D -> 3D (phi grid on the command line), map AMR data to level 1
  python3 dat_lift.py output/Star8_0034.dat --ndim 3 --level 1 \\
          --phi-min 0 --phi-max 0.1 --nx3 40 --bnx3 20 --periodic-phi
  # 1D -> 3D directly
  python3 dat_lift.py n1280_0100.dat --ndim 3 --par ../m1_Ostar_3D/parfiles/base.par
"""
import argparse
import math
import os
import re
import struct
import sys

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', '..', '..', 'tools', 'python'))
from amrvac_pytools.datfiles.reading.datfile_utilities import (  # noqa: E402
    get_header, get_tree_info)

SIZE_INT = 4
SIZE_DOUBLE = 8
NAME_LEN = 16
ZERO_IF_ABSENT = ('m2', 'm3', 'F2', 'F3')
DIM_NAME = {1: 'r', 2: 'theta', 3: 'phi'}


def default_vars(ndir):
    return (['rho'] + ['m%d' % i for i in range(1, ndir + 1)] + ['e', 'Er']
            + ['F%d' % i for i in range(1, ndir + 1)])


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def fortran_real(x):
    s = repr(float(x))
    return s.replace('e', 'd') if 'e' in s else s + 'd0'


def morton_encode(ig):
    """mortonEncode of mod_space_filling_curve.t for 0-based indices ig."""
    nd = len(ig)
    ans = 0
    for b in range(64 // nd + 1):
        for d, i in enumerate(ig):
            ans |= ((i >> b) & 1) << (nd * b + d)
    return ans


def level1_order(ng):
    """Level-1 blocks (1-based index tuples) in amrvac's Morton order: the
    in-domain blocks of the covering power-of-2 cube, ordered by Morton code."""
    grids = np.meshgrid(*[np.arange(1, n + 1) for n in ng], indexing='ij')
    blocks = [tuple(int(g) for g in t) for t in zip(*[g.ravel() for g in grids])]
    return sorted(blocks, key=lambda ig: morton_encode([i - 1 for i in ig]))


def read_par_values(files, keys):
    """Minimal Fortran-namelist reader: last value of each key (lower case)."""
    out = {}
    pat = re.compile(r'([A-Za-z_][\w()]*)\s*=\s*([^,=]+?)'
                     r'(?=\s*(?:,|$|\s[A-Za-z_][\w()]*\s*=))')
    for fn in files:
        with open(fn) as f:
            for line in f:
                for m in pat.finditer(line.split('!')[0].strip()):
                    k = m.group(1).lower()
                    if k in keys:
                        out[k] = m.group(2).strip()
    return out


def par_periodic(files, idim):
    """True/False if typeboundary_min/max<idim> in the par files says
    'periodic' / something else, None if not given (last occurrence wins)."""
    res = None
    pat = re.compile(r'typeboundary_(?:min|max)%d\s*=\s*(.*)' % idim, re.I)
    for fn in files:
        with open(fn) as f:
            for line in f:
                m = pat.search(line.split('!')[0])
                if m:
                    res = bool(re.search(r"['\"]periodic['\"]", m.group(1), re.I))
    return res


def par_float(s):
    return float(s.lower().replace('d', 'e'))


def par_int(s):
    return int(par_float(s))


def par_units_for(rad, hint=None):
    """A par-file value x with x*2.0*pi == rad exactly (amrvac computes
    xprob*two*dpi). Tries the hint, short decimals and neighbouring doubles."""
    cands = []
    if hint is not None:
        cands.append(hint)
    x0 = rad / (2.0 * math.pi)
    cands += [round(x0, k) for k in range(4, 18)] + [x0]
    x = x0
    for _ in range(4):
        x = np.nextafter(x, np.inf)
        cands.append(float(x))
    x = x0
    for _ in range(4):
        x = np.nextafter(x, -np.inf)
        cands.append(float(x))
    for c in cands:
        if c * 2.0 * math.pi == rad:
            return c
    return None


# --------------------------------------------------------------------------
# source snapshot
# --------------------------------------------------------------------------
def read_snapshot(fn):
    """Header, raw header bytes and {(lvl, ig tuple): w} with ghost cells
    stripped, w of shape block_nx + (nw,)."""
    with open(fn, 'rb') as f:
        h = get_header(f)
        nd = h['ndim']
        if h['datfile_version'] != 5:
            sys.exit('Only version-5 snapshots are supported (found %d)' % h['datfile_version'])
        if h.get('staggered', False):
            sys.exit('Staggered-grid snapshots are not supported')
        lvls, ixs, offs = get_tree_info(f)
        offs = offs - 2 * nd * SIZE_INT          # undo ghost-record skip
        f.seek(0)
        raw = f.read(h['offset_tree'])
        bnx = [int(b) for b in h['block_nx']]
        nw = h['nw']
        blocks = {}
        for lvl, ix, off in zip(lvls, ixs, offs):
            f.seek(off)
            g = np.frombuffer(f.read(2 * nd * SIZE_INT), dtype='<i4')
            lo, hi = g[:nd], g[nd:]
            shape = tuple(bnx[d] + lo[d] + hi[d] for d in range(nd)) + (nw,)
            w = np.frombuffer(f.read(int(np.prod(shape)) * SIZE_DOUBLE), dtype='<f8')
            w = w.reshape(shape, order='F')
            w = w[tuple(slice(lo[d], lo[d] + bnx[d]) for d in range(nd))]
            blocks[(int(lvl), tuple(int(i) for i in np.atleast_1d(ix)))] = w
    return h, raw, blocks


def cell_weights(h, lvl, nd, uniform_r):
    """Per-dimension volume weights of all cells at level lvl (spherical)."""
    n = [int(h['domain_nx'][d]) * 2 ** (lvl - 1) for d in range(nd)]
    ws = []
    for d in range(nd):
        lo, hi = h['xmin'][d], h['xmax'][d]
        if d == 0:
            e = (np.linspace(lo, hi, n[0] + 1) if uniform_r else
                 lo * ((hi / lo) ** (1.0 / n[0])) ** np.arange(n[0] + 1))
            ws.append(e[1:] ** 3 - e[:-1] ** 3)
        elif d == 1:
            e = np.linspace(lo, hi, n[1] + 1)
            ws.append(np.cos(e[:-1]) - np.cos(e[1:]))
        else:
            ws.append(np.full(n[d], (hi - lo) / n[d]))
    return ws


def map_to_level(h, blocks, L, uniform_r):
    """Uniform array (n1[, n2], nw) at level L."""
    nd = h['ndim']
    bnx = np.array([int(b) for b in h['block_nx']])
    nL = [int(h['domain_nx'][d]) * 2 ** (L - 1) for d in range(nd)]
    U = np.full(tuple(nL) + (h['nw'],), np.nan)
    wcache = {}
    for (lvl, ig), w in blocks.items():
        if lvl <= L:
            r = 2 ** (L - lvl)
            for d in range(nd):
                w = np.repeat(w, r, axis=d)
            j0 = [(ig[d] - 1) * bnx[d] * r for d in range(nd)]
            U[tuple(slice(j0[d], j0[d] + bnx[d] * r) for d in range(nd))] = w
        else:
            r = 2 ** (lvl - L)
            if lvl not in wcache:
                wcache[lvl] = cell_weights(h, lvl, nd, uniform_r)
            wd = [wcache[lvl][d][(ig[d] - 1) * bnx[d]:ig[d] * bnx[d]] for d in range(nd)]
            V = wd[0]
            for d in range(1, nd):
                V = np.multiply.outer(V, wd[d])
            shape = []
            for d in range(nd):
                shape += [bnx[d] // r, r]
            ax = tuple(2 * d + 1 for d in range(nd))
            num = (w * V[..., None]).reshape(shape + [w.shape[-1]]).sum(axis=ax)
            den = V.reshape(shape).sum(axis=ax)[..., None]
            j0 = [((ig[d] - 1) * bnx[d]) // r for d in range(nd)]
            U[tuple(slice(j0[d], j0[d] + bnx[d] // r) for d in range(nd))] = num / den
    if np.isnan(U).any():
        sys.exit('Could not fill the source data at level %d (inconsistent tree?)' % L)
    return U


# --------------------------------------------------------------------------
# target snapshot
# --------------------------------------------------------------------------
def name_bytes(s):
    b = s.encode()
    if len(b) > NAME_LEN:
        sys.exit('Name too long for a .dat header: %s' % s)
    return b.ljust(NAME_LEN)


def source_header_parts(h, raw):
    """Periodic flags, geometry and stagger bytes and the physics tail of the
    source header (layout of snapshot_write_header1)."""
    nd, nw = h['ndim'], h['nw']
    pos_per = 48 + 16 * nd + 8 * nd
    periodic = list(struct.unpack_from('<%di' % nd, raw, pos_per))
    pos_geom = pos_per + 4 * nd
    geom = raw[pos_geom:pos_geom + NAME_LEN]
    stagger = raw[pos_geom + NAME_LEN:pos_geom + NAME_LEN + 4]
    tail = bytearray(raw[pos_geom + NAME_LEN + 4 + NAME_LEN * nw:])
    return periodic, geom, stagger, tail


def build_header(h, raw, nd, ndir, names, it, time, xmin, xmax, domain_nx,
                 block_nx, periodic, nleafs, offset_tree, offset_blocks, reset):
    _, geom, stagger, tail = source_header_parts(h, raw)
    if reset:
        # snapshotnext, slicenext, collapsenext. A real snapshot 0000 stores
        # snapshotnext=1; on a restart at it=0 amrvac subtracts one again
        # (amrvac.t) and re-saves the start as 0000
        struct.pack_into('<iii', tail, len(tail) - 12, 1, 0, 0)
    hdr = bytearray()
    hdr += struct.pack('<10i', h['datfile_version'], offset_tree, offset_blocks,
                       len(names), ndir, nd, 1, nleafs, 0, it)
    hdr += struct.pack('<d', time)
    hdr += struct.pack('<%dd' % nd, *xmin) + struct.pack('<%dd' % nd, *xmax)
    hdr += struct.pack('<%di' % nd, *domain_nx) + struct.pack('<%di' % nd, *block_nx)
    hdr += struct.pack('<%di' % nd, *[int(bool(p)) for p in periodic])
    hdr += geom + stagger
    for n in names:
        hdr += name_bytes(n)
    hdr += tail
    return bytes(hdr)


def main(argv=None, default_ndim=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('datfile', help='source snapshot (.dat), 1D or 2D spherical')
    ap.add_argument('--ndim', type=int, choices=(2, 3), default=default_ndim,
                    required=default_ndim is None, help='target dimensionality')
    ap.add_argument('--ndir', type=int, choices=(2, 3),
                    help='target ndir (default: target ndim; 3 for 2.5D)')
    ap.add_argument('--par', nargs='+', default=[],
                    help='target par file(s) for the grid of the new dimensions '
                    '(and block sizes)')
    ap.add_argument('--theta-min', type=float, help='xprobmin2 (par units)')
    ap.add_argument('--theta-max', type=float, help='xprobmax2 (par units)')
    ap.add_argument('--nx2', type=int, help='domain_nx2 (new theta only)')
    ap.add_argument('--bnx2', type=int, help='block_nx2')
    ap.add_argument('--phi-min', type=float, help='xprobmin3 (par units)')
    ap.add_argument('--phi-max', type=float, help='xprobmax3 (par units)')
    ap.add_argument('--nx3', type=int, help='domain_nx3')
    ap.add_argument('--bnx3', type=int, help='block_nx3')
    g = ap.add_mutually_exclusive_group()
    g.add_argument('--periodic-phi', dest='periodic_phi', action='store_true',
                   default=None, help='phi is periodic in the target run')
    g.add_argument('--no-periodic-phi', dest='periodic_phi', action='store_false',
                   help='phi is not periodic in the target run')
    ap.add_argument('--bnx1', type=int, help='block_nx1 of the target run')
    ap.add_argument('--level', type=int, help='level to map source AMR data to '
                    '(default: finest level in the source)')
    ap.add_argument('--uniform-r', action='store_true',
                    help='source radial grid is uniform, not stretched (only '
                    'used to volume-average finer AMR levels)')
    ap.add_argument('--vars', nargs='+', help='output variables, in target order')
    ap.add_argument('--reset', action='store_true',
                    help='start the target run at it=0, t=0, numbering from 0000')
    ap.add_argument('-o', '--output', help='output .dat (a 4-digit index is '
                    'appended if missing, as amrvac requires)')
    ap.add_argument('--out-par', help='output par file (default: <output>.par)')
    args = ap.parse_args(argv)

    h, raw, blocks = read_snapshot(args.datfile)
    s, t = h['ndim'], args.ndim
    if not h['geometry'].lower().startswith('spherical'):
        sys.exit('Only spherical snapshots are supported (geometry=%s)' % h['geometry'])
    if t <= s:
        sys.exit('Target ndim (%d) must be larger than the source ndim (%d)' % (t, s))
    ndir = args.ndir or t
    if ndir < t:
        sys.exit('--ndir must be at least the target ndim')
    names_out = args.vars or default_vars(ndir)

    # ---- par files -----------------------------------------------------------
    keys = ['xprobmin%d' % d for d in (1, 2, 3)] + ['xprobmax%d' % d for d in (1, 2, 3)] \
        + ['domain_nx%d' % d for d in (1, 2, 3)] + ['block_nx%d' % d for d in (1, 2, 3)] \
        + ['qstretch_baselevel', 'qstretch_baselevel(1)']
    pv = read_par_values(args.par, keys) if args.par else {}

    def pick(cli, key, conv):
        if cli is not None:
            return cli
        return conv(pv[key]) if key in pv else None

    cli_min = {2: args.theta_min, 3: args.phi_min}
    cli_max = {2: args.theta_max, 3: args.phi_max}
    cli_nx = {2: args.nx2, 3: args.nx3}
    cli_bnx = {1: args.bnx1, 2: args.bnx2, 3: args.bnx3}

    # ---- source mapping (kept dimensions) --------------------------------------
    L = args.level or h['levmax']
    if L > h['levmax'] and h['levmax'] == 1:
        pass   # refining a single-level source by copying is allowed
    warnings = []
    xmin, xmax, xmin_par, xmax_par, nx, bnx, periodic, origin = [], [], [], [], [], [], [], []
    src_periodic = source_header_parts(h, raw)[0]
    for d in range(1, s + 1):
        lo, hi = float(h['xmin'][d - 1]), float(h['xmax'][d - 1])
        n = int(h['domain_nx'][d - 1]) * 2 ** (L - 1)
        if d == 1:
            lo_p, hi_p = lo, hi
            for k, v in (('xprobmin1', lo), ('xprobmax1', hi)):
                if k in pv and par_float(pv[k]) != v:
                    warnings.append('%s=%s in par file differs from the source (%.17g)'
                                    % (k, pv[k], v))
        else:
            hint_lo = par_float(pv['xprobmin%d' % d]) if 'xprobmin%d' % d in pv else None
            hint_hi = par_float(pv['xprobmax%d' % d]) if 'xprobmax%d' % d in pv else None
            lo_p, hi_p = par_units_for(lo, hint_lo), par_units_for(hi, hint_hi)
            if lo_p is None or hi_p is None:
                sys.exit('No par-file value reproduces the source %s range exactly '
                         '(x*2*pi); please report' % DIM_NAME[d])
            for k, hint, v in (('xprobmin%d' % d, hint_lo, lo_p), ('xprobmax%d' % d, hint_hi, hi_p)):
                if hint is not None and hint * 2.0 * math.pi != (lo if 'min' in k else hi):
                    warnings.append('%s=%s in par file differs from the source (%.17g par units)'
                                    % (k, pv[k], v))
        if 'domain_nx%d' % d in pv and par_int(pv['domain_nx%d' % d]) != n:
            warnings.append('domain_nx%d=%s in par file; the output has %d (the generated '
                            'par file sets it)' % (d, pv['domain_nx%d' % d], n))
        b = pick(cli_bnx[d], 'block_nx%d' % d, par_int) or int(h['block_nx'][d - 1])
        xmin.append(lo); xmax.append(hi); xmin_par.append(lo_p); xmax_par.append(hi_p)
        nx.append(n); bnx.append(b); periodic.append(bool(src_periodic[d - 1]))
        origin.append('source')
    if any(k in pv for k in ('qstretch_baselevel', 'qstretch_baselevel(1)')):
        warnings.append('qstretch_baselevel is set in the par file: make sure the source '
                        'run used the same radial stretching')
    if h['levmax'] > 1:
        warnings.append('source has AMR (levmax=%d), mapped to level %d: kept dimensions '
                        'have %s cells' % (h['levmax'], L, ' x '.join(map(str, nx))))

    # ---- new dimensions --------------------------------------------------------
    for d in range(s + 1, t + 1):
        a_par = pick(cli_min[d], 'xprobmin%d' % d, par_float)
        b_par = pick(cli_max[d], 'xprobmax%d' % d, par_float)
        n = pick(cli_nx[d], 'domain_nx%d' % d, par_int)
        b = pick(cli_bnx[d], 'block_nx%d' % d, par_int)
        miss = [k for k, v in (('xprobmin%d' % d, a_par), ('xprobmax%d' % d, b_par),
                               ('domain_nx%d' % d, n), ('block_nx%d' % d, b)) if v is None]
        if miss:
            sys.exit('%s grid incomplete, missing: %s (give --par and/or the command-line '
                     'options)' % (DIM_NAME[d], ', '.join(miss)))
        if not b_par > a_par:
            sys.exit('xprobmax%d must be larger than xprobmin%d' % (d, d))
        # amrvac multiplies xprob^LIM2,3 by two*dpi for spherical geometry;
        # same operation order so the header matches the par file bit for bit
        xmin.append(a_par * 2.0 * math.pi); xmax.append(b_par * 2.0 * math.pi)
        xmin_par.append(a_par); xmax_par.append(b_par)
        nx.append(n); bnx.append(b); origin.append('new')
        if d == 3 and args.periodic_phi is not None:
            per = args.periodic_phi
        else:
            per = par_periodic(args.par, d)
            if per is None:
                per = False
                warnings.append('periodicity of %s not given (--periodic-phi or '
                                'typeboundary_min%d in --par); written as non-periodic. '
                                'It must match the target run.' % (DIM_NAME[d], d))
        periodic.append(per)

    for d in range(t):
        if nx[d] % bnx[d]:
            sys.exit('domain_nx%d=%d is not a multiple of block_nx%d=%d'
                     % (d + 1, nx[d], d + 1, bnx[d]))

    # ---- variables -------------------------------------------------------------
    U = map_to_level(h, blocks, L, args.uniform_r)
    del blocks
    cols = []
    for v in names_out:
        if v in h['w_names']:
            cols.append(U[..., h['w_names'].index(v)])
        elif v in ZERO_IF_ABSENT:
            cols.append(np.zeros(U.shape[:-1]))
        else:
            sys.exit('Variable %s not in the source (%s) and not a zero-filled vector '
                     'component' % (v, ' '.join(h['w_names'])))
    P = np.stack(cols, axis=-1)          # shape: kept-dims + (nw_out,)
    del U
    nw_out = P.shape[-1]

    # ---- output names ----------------------------------------------------------
    m = re.match(r'^(.*?)(\d{4})\.dat$', os.path.basename(args.datfile))
    idx = '0000' if (args.reset or not m) else m.group(2)
    if args.output:
        out = args.output
        if not re.search(r'\d{4}\.dat$', out):
            out = os.path.splitext(out)[0] + idx + '.dat'
    else:
        base = m.group(1) if m else os.path.splitext(os.path.basename(args.datfile))[0]
        out = os.path.join(os.path.dirname(args.datfile), '%s%dd%s.dat' % (base, t, idx))
    parfile = args.out_par or os.path.splitext(out)[0] + '.par'

    # ---- write -----------------------------------------------------------------
    ng = [nx[d] // bnx[d] for d in range(t)]
    leaves = level1_order(ng)
    nleafs = len(leaves)
    it = 0 if args.reset else h['it']
    time = 0.0 if args.reset else h['time']

    def header(ot, ob):
        return build_header(h, raw, t, ndir, names_out, it, time, xmin, xmax, nx, bnx,
                            periodic, nleafs, ot, ob, args.reset)

    offset_tree = len(header(0, 0))
    offset_blocks = offset_tree + SIZE_INT * nleafs * (2 + t) + 8 * nleafs
    hdr = header(offset_tree, offset_blocks)
    block_bytes = 2 * t * SIZE_INT + int(np.prod(bnx)) * nw_out * SIZE_DOUBLE
    offsets = offset_blocks + block_bytes * np.arange(nleafs, dtype=np.int64)

    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    full_shape = tuple(bnx) + (nw_out,)
    with open(out, 'wb') as f:
        f.write(hdr)
        f.write(np.ones(nleafs, dtype='<i4').tobytes())          # leaf flags
        f.write(np.ones(nleafs, dtype='<i4').tobytes())          # levels
        f.write(np.array(leaves, dtype='<i4').tobytes())         # block indices
        f.write(offsets.astype('<i8').tobytes())
        assert f.tell() == offset_blocks
        for ig in leaves:
            sub = P[tuple(slice((ig[d] - 1) * bnx[d], ig[d] * bnx[d]) for d in range(s))]
            sub = sub.reshape(sub.shape[:s] + (1,) * (t - s) + (nw_out,))
            w = np.broadcast_to(sub, full_shape)
            f.write(np.zeros(2 * t, dtype='<i4').tobytes())      # no ghost cells
            f.write(np.asarray(w, dtype='<f8').tobytes(order='F'))
        assert f.tell() == offset_blocks + block_bytes * nleafs

    with open(parfile, 'w') as f:
        f.write("! Generated by dat_lift.py from %s (%dD -> %dD)\n" % (args.datfile, s, t))
        f.write("! Give this file last: ./amrvac -i <base>.par <star>.par %s\n" % parfile)
        f.write("! - radial stretching (stretch_dim(1), qstretch_baselevel) must be the\n"
                "!   same as in the source run\n")
        f.write("! - periodicity written: %s; the boundary types of the target run must\n"
                "!   match (typeboundary_min/max<d> = 'periodic' for periodic dimensions)\n"
                % ', '.join('%s %s' % (DIM_NAME[d + 1], 'periodic' if periodic[d] else
                                       'non-periodic') for d in range(t)))
        f.write("! - ndir=%d, variables %s\n\n" % (ndir, ' '.join(names_out)))
        f.write("&filelist\n  restart_from_file = '%s'\n/\n\n&meshlist\n" % out)
        for d in range(t):
            f.write("  domain_nx%d = %d\n  block_nx%d  = %d\n  xprobmin%d  = %s\n"
                    "  xprobmax%d  = %s\n" % (d + 1, nx[d], d + 1, bnx[d], d + 1,
                                              fortran_real(xmin_par[d]), d + 1,
                                              fortran_real(xmax_par[d])))
        f.write("/\n")

    print('Input : %s  %dD, domain_nx=%s, levmax=%d, it=%d, t=%.6g'
          % (args.datfile, s, 'x'.join(str(int(x)) for x in h['domain_nx']),
             h['levmax'], h['it'], h['time']))
    print('        vars %s' % ' '.join(h['w_names']))
    print('Output: %s  %dD, ndir=%d, %s cells, blocks %s, %d leaves at level 1, it=%d, t=%.6g'
          % (out, t, ndir, ' x '.join(map(str, nx)), ' x '.join(map(str, bnx)),
             nleafs, it, time))
    for d in range(t):
        print('        %-5s [%.10g, %.10g] (par units), %s, %s'
              % (DIM_NAME[d + 1], xmin_par[d], xmax_par[d],
                 'periodic' if periodic[d] else 'non-periodic', origin[d]))
    print('        vars %s' % ' '.join(names_out))
    print('Par file: %s' % parfile)
    for w_ in warnings:
        print('WARNING: ' + w_)


if __name__ == '__main__':
    main()
