#!/usr/bin/env python3
"""
Expand the theta range of a 2D (r, theta) MPI-AMRVAC snapshot by mirroring.

The data of the original wedge [t0, t1] is copied into the new range [a, b]
by repeated reflection about the original theta boundaries: e.g. for an
original wedge 0.2...0.3 and a new range 0.0...0.5, theta in [0.0, 0.2] is the
mirror image of [0.2, 0.3] about 0.2, and [0.3, 0.5] the mirror image about
0.3 (and so on for wider extensions). Mirrored copies have their theta-vector
components (m2, F2, ...) negated, consistent with 'asymm' theta boundaries.

Two modes:
  * default (exact): whole AMR leaf blocks are mirrored, so the refinement
    structure is preserved exactly and no interpolation is done. The extension
    widths (t0-a) and (b-t1) must be integer multiples of the level-1 block
    width in theta.
  * --resample: for arbitrary new ranges. The snapshot is first mapped onto a
    uniform grid at one refinement level (--level, default: the file's levmax;
    coarser blocks are copied piecewise-constant, finer ones averaged), the
    mirrored theta profile is then linearly interpolated onto the new theta
    cells (domain_nx2 from --nx2, default: the multiple of block_nx2 closest to
    the original theta resolution), and all output blocks are at that level.
    Values are interpolated in theta only; r is untouched.

Outputs:
  * a new .dat snapshot (version-5 format, as written by write_snapshot)
  * a small .par file, to be given LAST on the command line, e.g.
        ./amrvac -i parfiles/Ostar.par restart_expand.par

Usage:
  python expand_theta.py output/Ostar0010.dat 0.0 0.5 \
         [-o out.dat] [--par restart_expand.par] [--flip m2 F2]
         [--resample [--nx2 N] [--level L]]
"""
import argparse
import math
import re
import os
import struct
import sys

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', '..', '..', 'tools', 'python'))
from amrvac_pytools.datfiles.reading.datfile_utilities import (  # noqa: E402
    get_header, get_tree_info)

SIZE_INT = 4
SIZE_DOUBLE = 8
# theta-odd quantities: theta components of vectors, r-theta tensor
# components, and (for this setup) the angular momentum am = r*m2
DEFAULT_FLIP = ('m2', 'v2', 'b2', 'F2', 'Fr2', 'f2', 'P_rt', 'am')


def morton_encode_2d(i, j):
    """Interleave bits as mortonEncode in mod_space_filling_curve.t."""
    ans = 0
    for b in range(32):
        ans |= ((i >> b) & 1) << (2 * b)
        ans |= ((j >> b) & 1) << (2 * b + 1)
    return ans


def level1_order(ng1, ng2):
    """Level-1 blocks (1-based) in amrvac's Morton order.

    amrvac numbers the blocks of the covering power-of-2 square with a Morton
    curve and removes the out-of-domain ones, i.e. the in-domain blocks are
    ordered by their Morton code."""
    blocks = [(i1, i2) for i2 in range(1, ng2 + 1) for i1 in range(1, ng1 + 1)]
    return sorted(blocks, key=lambda ig: morton_encode_2d(ig[0] - 1, ig[1] - 1))


def fortran_real(x):
    s = repr(float(x))
    return s.replace('e', 'd') if 'e' in s else s + 'd0'


def read_snapshot(fn):
    """Return header dict, raw header bytes and {(lvl, i1, i2): (w, g1lo, g1hi)}.

    theta ghost cells are stripped from w; r ghost counts are returned."""
    with open(fn, 'rb') as f:
        h = get_header(f)
        if h['ndim'] != 2:
            sys.exit('Only 2D snapshots are supported (ndim=%d)' % h['ndim'])
        if h.get('staggered', False):
            sys.exit('Staggered-grid snapshots are not supported')
        lvls, ixs, offs = get_tree_info(f)
        offs = offs - 2 * h['ndim'] * SIZE_INT  # undo ghost-record skip
        f.seek(0)
        raw_header = bytearray(f.read(h['offset_tree']))

        nw = h['nw']
        bnx = h['block_nx']
        blocks = {}
        for lvl, ix, off in zip(lvls, ixs, offs):
            f.seek(off)
            ng = np.frombuffer(f.read(4 * SIZE_INT), dtype='<i4')
            shape = (bnx[0] + ng[0] + ng[2], bnx[1] + ng[1] + ng[3], nw)
            n = int(np.prod(shape))
            w = np.frombuffer(f.read(n * SIZE_DOUBLE), dtype='<f8')
            w = w.reshape(shape, order='F')
            w = w[:, ng[1]:ng[1] + bnx[1], :]
            blocks[(int(lvl), int(ix[0]), int(ix[1]))] = (w, int(ng[0]), int(ng[2]))
    return h, raw_header, blocks


def build_forest(ng1, ng2, is_leaf, max_level):
    """Traverse the forest in amrvac order (write_node: ic2 outer, ic1 inner).

    is_leaf(lvl, i1, i2) decides whether a node is a leaf.
    Returns leaf flags (per node), leaves [(lvl, i1, i2)], nparents."""
    flags, leaves = [], []
    nparents = 0

    def visit(lvl, i1, i2):
        nonlocal nparents
        if is_leaf(lvl, i1, i2):
            flags.append(1)
            leaves.append((lvl, i1, i2))
        else:
            if lvl >= max_level:
                raise RuntimeError('inconsistent forest at (%d,%d,%d)' % (lvl, i1, i2))
            flags.append(0)
            nparents += 1
            for ic2 in (1, 2):
                for ic1 in (1, 2):
                    visit(lvl + 1, 2 * (i1 - 1) + ic1, 2 * (i2 - 1) + ic2)

    for i1, i2 in level1_order(ng1, ng2):
        visit(1, i1, i2)
    return flags, leaves, nparents


def write_snapshot(out, raw_header, flags, leaves, get_block,
                   xmin2, xmax2, domain_nx2):
    """Write a snapshot. get_block(leaf) returns (w, g1lo, g1hi) with
    w of shape (nx1+g1lo+g1hi, nx2, nw) (no theta ghost cells)."""
    raw_header = bytearray(raw_header)
    nleafs = len(leaves)
    nparents = len(flags) - nleafs
    levmax = max(l[0] for l in leaves)
    offset_tree = len(raw_header)
    offset_blocks = (offset_tree + SIZE_INT * len(flags)
                     + SIZE_INT * nleafs * 3 + 8 * nleafs)

    # Patch header fields (layout of snapshot_write_header1)
    struct.pack_into('<ii', raw_header, 4, offset_tree, offset_blocks)
    struct.pack_into('<iii', raw_header, 24, levmax, nleafs, nparents)
    struct.pack_into('<d', raw_header, 48 + 8, xmin2)
    struct.pack_into('<d', raw_header, 64 + 8, xmax2)
    struct.pack_into('<i', raw_header, 80 + 4, domain_nx2)

    # Block offsets (record = 4 ghost counts + data)
    offsets = []
    pos = offset_blocks
    for leaf in leaves:
        w = get_block(leaf)[0]
        offsets.append(pos)
        pos += 4 * SIZE_INT + w.size * SIZE_DOUBLE

    with open(out, 'wb') as f:
        f.write(raw_header)
        f.write(np.array(flags, dtype='<i4').tobytes())
        f.write(np.array([l[0] for l in leaves], dtype='<i4').tobytes())
        f.write(np.array([[l[1], l[2]] for l in leaves], dtype='<i4').tobytes())
        f.write(np.array(offsets, dtype='<i8').tobytes())
        assert f.tell() == offset_blocks
        for leaf in leaves:
            w, g1lo, g1hi = get_block(leaf)
            f.write(np.array([g1lo, 0, g1hi, 0], dtype='<i4').tobytes())
            f.write(np.asarray(w, dtype='<f8').tobytes(order='F'))
        assert f.tell() == pos
    return levmax, nleafs


def uniform_level_data(h, blocks, L):
    """Map all blocks onto a uniform array at refinement level L (no ghosts)."""
    bnx = h['block_nx']
    nx = h['domain_nx'] * 2 ** (L - 1)
    U = np.empty((nx[0], nx[1], h['nw']))
    for (lvl, i1, i2), (w, g1lo, g1hi) in blocks.items():
        w = w[g1lo:g1lo + bnx[0]]
        if lvl <= L:
            r = 2 ** (L - lvl)
            w = w.repeat(r, axis=0).repeat(r, axis=1)
            n = bnx * r
            j0 = np.array([i1 - 1, i2 - 1]) * n
        else:
            r = 2 ** (lvl - L)
            w = w.reshape(bnx[0] // r, r, bnx[1] // r, r, -1).mean(axis=(1, 3))
            n = bnx // r
            j0 = (np.array([i1 - 1, i2 - 1]) // r) * bnx + \
                (np.array([i1 - 1, i2 - 1]) % r) * n
        U[j0[0]:j0[0] + n[0], j0[1]:j0[1] + n[1]] = w
    return U


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('datfile', help='input .dat snapshot')
    ap.add_argument('theta_min', type=float,
                    help='new xprobmin2, in par-file units (see --radians)')
    ap.add_argument('theta_max', type=float,
                    help='new xprobmax2, in par-file units (see --radians)')
    ap.add_argument('--radians', action='store_true',
                    help='theta_min/max are given in radians. By default they '
                    'are in par-file units: for spherical geometry amrvac '
                    'multiplies xprobmin2/xprobmax2 by 2*pi')
    ap.add_argument('-o', '--output', help='output .dat file '
                    '(default: <base>theta<NNNN>.dat; a 4-digit index is '
                    'appended if missing, as amrvac requires)')
    ap.add_argument('--par', help='output par file (default: <output>.par)')
    ap.add_argument('--flip', nargs='*', default=None,
                    help='variables negated in mirrored copies '
                    '(default: those among %s)' % ' '.join(DEFAULT_FLIP))
    ap.add_argument('--resample', action='store_true',
                    help='allow arbitrary new ranges by interpolating in theta '
                    'onto a uniform grid at one level (see above)')
    ap.add_argument('--nx2', type=int, default=None,
                    help='with --resample: new domain_nx2 (multiple of block_nx2)')
    ap.add_argument('--level', type=int, default=None,
                    help='with --resample: refinement level of the output '
                    '(default: levmax of the input)')
    ap.add_argument('--tol', type=float, default=1e-8,
                    help='relative tolerance for block alignment check')
    args = ap.parse_args()

    # amrvac parses the snapshot index from restart_from_file
    # (get_snapshot_index: the 4 characters before '.dat'), so keep that form
    m = re.match(r'^(.*?)(\d{4})\.dat$', args.datfile)
    idx = m.group(2) if m else '0000'
    if args.output:
        out = args.output
        if not re.search(r'\d{4}\.dat$', out):
            out = os.path.splitext(out)[0] + idx + '.dat'
    else:
        base = m.group(1) if m else os.path.splitext(args.datfile)[0]
        out = base + 'theta' + idx + '.dat'
    parfile = args.par or os.path.splitext(out)[0] + '.par'
    a, b = args.theta_min, args.theta_max

    h, raw_header, blocks = read_snapshot(args.datfile)

    # amrvac (read_par_files) multiplies xprob^LIM2 by 2*pi for spherical
    # geometry; the snapshot header stores the scaled values (radians)
    scale = 2.0 * math.pi if h['geometry'].lower().startswith('spherical') else 1.0
    if args.radians:
        a_par, b_par = a / scale, b / scale
    else:
        a_par, b_par = a, b
        if scale != 1.0:   # same operation order as Fortran: x*two*dpi
            a, b = a * 2.0 * math.pi, b * 2.0 * math.pi
    if not b > a:
        sys.exit('theta_max must be larger than theta_min')

    t0, t1 = h['xmin'][1], h['xmax'][1]
    W = t1 - t0
    bnx = h['block_nx']
    dnx = h['domain_nx']
    ng1 = dnx[0] // bnx[0]
    ng2_old = dnx[1] // bnx[1]
    dblk = bnx[1] * W / dnx[1]          # level-1 block width in theta

    names = h['w_names']
    flip = [v for v in names if v in DEFAULT_FLIP] if args.flip is None else args.flip
    unknown = [v for v in flip if v not in names]
    if unknown:
        sys.exit('Unknown variable(s) in --flip: %s (available: %s)'
                 % (unknown, names))
    iflip = [names.index(v) for v in flip]

    if not args.resample:
        klo = (t0 - a) / dblk
        khi = (b - t1) / dblk
        tol = args.tol * max(1.0, abs(klo), abs(khi))
        if (klo < -tol or khi < -tol or abs(klo - round(klo)) > tol
                or abs(khi - round(khi)) > tol):
            sys.exit('New range [%.10g, %.10g] must extend the old range '
                     '[%.10g, %.10g] by integer multiples of the level-1 block '
                     'width %.10g (par-file units; x2*pi = radians for '
                     'spherical).\nNearest valid range: [%.10g, %.10g]\n'
                     'Or use --resample to interpolate onto the requested range.'
                     % (a_par, b_par, t0 / scale, t1 / scale, dblk / scale,
                        (t0 - max(0, round(klo)) * dblk) / scale,
                        (t1 + max(0, round(khi)) * dblk) / scale))
        klo, khi = int(round(klo)), int(round(khi))
        ng2_new = ng2_old + klo + khi
        domain_nx2_new = ng2_new * bnx[1]

        def source(lvl, j):
            """New 1-based theta block index at level lvl -> (old index, flipped)."""
            N = ng2_old * 2 ** (lvl - 1)
            q = (j - 1 - klo * 2 ** (lvl - 1)) % (2 * N)
            return (q + 1, False) if q < N else (2 * N - q, True)

        def is_leaf(lvl, i1, i2):
            return (lvl, i1, source(lvl, i2)[0]) in blocks

        def get_block(leaf):
            lvl, i1, i2 = leaf
            src, flipped = source(lvl, i2)
            w, g1lo, g1hi = blocks[(lvl, i1, src)]
            if flipped:
                w = w[:, ::-1, :].copy()
                w[:, :, iflip] *= -1.0
            return w, g1lo, g1hi

        flags, leaves, _ = build_forest(ng1, ng2_new, is_leaf, h['levmax'])
        note = '+%d / +%d level-1 blocks below/above, exact copy' % (klo, khi)

    else:
        L = args.level or h['levmax']
        dth1_old = W / dnx[1]
        if args.nx2 is None:
            nb = max(1, int(round((b - a) / dth1_old / bnx[1])))
            domain_nx2_new = nb * bnx[1]
        else:
            domain_nx2_new = args.nx2
        if domain_nx2_new % bnx[1] != 0:
            sys.exit('--nx2 must be a multiple of block_nx2=%d' % bnx[1])
        ng2_new = domain_nx2_new // bnx[1]

        U = uniform_level_data(h, blocks, L)
        del blocks
        n_old = U.shape[1]
        dth_old = W / n_old
        sgn = np.ones(h['nw'])
        sgn[iflip] = -1.0
        # one mirror ghost cell on each side (centres at t0 - dth/2, t1 + dth/2)
        Up = np.concatenate([U[:, :1] * sgn, U, U[:, -1:] * sgn], axis=1)
        del U

        n_new = domain_nx2_new * 2 ** (L - 1)
        th = a + (np.arange(n_new) + 0.5) * (b - a) / n_new
        u = np.mod(th - t0, 2 * W)
        mirrored = u > W
        u = np.where(mirrored, 2 * W - u, u)
        s = u / dth_old + 0.5                     # fractional index in Up
        k0 = np.clip(np.floor(s).astype(int), 0, n_old)
        fr = (s - k0)[None, :, None]
        Unew = (1.0 - fr) * Up[:, k0] + fr * Up[:, k0 + 1]
        del Up
        for i in iflip:
            Unew[:, mirrored, i] *= -1.0

        def is_leaf(lvl, i1, i2):
            return lvl == L

        def get_block(leaf):
            lvl, i1, i2 = leaf
            return (Unew[(i1 - 1) * bnx[0]:i1 * bnx[0],
                         (i2 - 1) * bnx[1]:i2 * bnx[1]], 0, 0)

        flags, leaves, _ = build_forest(ng1, ng2_new, is_leaf, L)
        note = ('resampled at level %d: dtheta %.6g -> %.6g (par units), '
                'interpolated in theta; amrvac will re-coarsen at its first regrid'
                % (L, dth_old / scale, (b - a) / n_new / scale))

    levmax, nleafs = write_snapshot(out, raw_header, flags, leaves, get_block,
                                    a, b, domain_nx2_new)

    with open(parfile, 'w') as f:
        f.write("! Generated by expand_theta.py from %s\n" % args.datfile)
        f.write("! Give this file last: ./amrvac -i <your>.par %s\n\n" % parfile)
        f.write("&filelist\n  restart_from_file = '%s'\n/\n\n" % out)
        f.write("&meshlist\n  xprobmin2  = %s\n  xprobmax2  = %s\n"
                "  domain_nx2 = %d\n/\n" % (fortran_real(a_par), fortran_real(b_par),
                                            domain_nx2_new))

    print('Input : %s  xprob2 [%.10g, %.10g] (theta [%.10g, %.10g] rad), '
          'domain_nx2=%d, %d leaves'
          % (args.datfile, t0 / scale, t1 / scale, t0, t1, dnx[1], h['nleafs']))
    print('Output: %s  xprob2 [%.10g, %.10g] (theta [%.10g, %.10g] rad), '
          'domain_nx2=%d, %d leaves' % (out, a_par, b_par, a, b,
                                         domain_nx2_new, nleafs))
    print('        (%s)' % note)
    print('Negated in mirrored copies: %s' % (', '.join(flip) or 'none'))
    print('Par file: %s' % parfile)
    print('Reminders: refine_max_level >= %d; it=%d, time=%.6g are kept '
          '(adjust it_max/time_max or use reset_it/reset_time).'
          % (levmax, h['it'], h['time']))


if __name__ == '__main__':
    main()
