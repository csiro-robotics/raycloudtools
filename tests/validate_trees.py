#!/usr/bin/env python3
"""
Validate rayextract trees output quality against the input point cloud.

Checks:
  - Correct tree count
  - All input points present in segmented output
  - Mesh height completeness (mesh reaches near the top of the cloud)
  - Coverage: fraction of input points close to the mesh surface
  - Spatial extent: mesh XY footprint covers the cloud footprint

Usage:
  python3 tests/validate_trees.py \\
      --binary build/bin/rayextract \\
      --cloud  test_data/split_trees/1.ply \\
      --mesh   test_data/split_trees/1_mesh.ply \\
      [--expected-trees 1] \\
      [--min-coverage 0.70] \\
      [--no-run]      # skip extraction, validate existing output files

Exit code 0 = all checks passed.
"""

import sys
import os
import struct
import re
import math
import subprocess
import argparse
from collections import defaultdict


# ── PLY helpers ──────────────────────────────────────────────────────────────

_TYPE = {
    'float': ('f', 4), 'double': ('d', 8),
    'int': ('i', 4), 'uint': ('I', 4),
    'short': ('h', 2), 'ushort': ('H', 2),
    'uchar': ('B', 1), 'char': ('b', 1),
}

def _ply_header(path):
    with open(path, 'rb') as f:
        hdr = b''
        while True:
            line = f.readline()
            hdr += line
            if line.strip() == b'end_header':
                break
    return hdr.decode('ascii', errors='replace')


def _vertex_props(hdr):
    """Return (nv, props) where props is the list of (type, name) for the vertex element only."""
    nv = 0
    props = []
    in_vertex = False
    for line in hdr.splitlines():
        line = line.strip()
        m = re.match(r'element (\w+) (\d+)', line)
        if m:
            in_vertex = (m.group(1) == 'vertex')
            if in_vertex:
                nv = int(m.group(2))
            continue
        if in_vertex and line.startswith('property '):
            parts = line.split()
            if len(parts) >= 3 and parts[1] != 'list':
                props.append((parts[1], parts[2]))
    return nv, props


def read_ply_xyz(path, sample_every=100):
    """Return sampled XYZ positions from a PLY vertex list."""
    hdr = _ply_header(path)
    nv, props = _vertex_props(hdr)

    fmt = '<' + ''.join(_TYPE.get(t, ('B', 1))[0] for t, _ in props)
    stride = struct.calcsize(fmt)
    xi = next(i for i, (_, n) in enumerate(props) if n == 'x')
    yi = next(i for i, (_, n) in enumerate(props) if n == 'y')
    zi = next(i for i, (_, n) in enumerate(props) if n == 'z')

    with open(path, 'rb') as f:
        while f.readline().strip() != b'end_header':
            pass
        points = []
        for i in range(nv):
            raw = f.read(stride)
            if len(raw) < stride:
                break
            if i % sample_every == 0:
                row = struct.unpack_from(fmt, raw)
                points.append((row[xi], row[yi], row[zi]))
    return points, nv


def read_ply_vertex_count(path):
    hdr = _ply_header(path)
    m = re.search(r'element vertex (\d+)', hdr)
    return int(m.group(1)) if m else 0


# ── Spatial hash ─────────────────────────────────────────────────────────────

class VoxelGrid:
    """Fast 3-D nearest-neighbour via voxel hashing."""
    def __init__(self, pts, cell_size=2.0):
        self.cell = cell_size
        self.grid = defaultdict(list)
        for i, (x, y, z) in enumerate(pts):
            self.grid[self._key(x, y, z)].append(i)
        self.pts = pts

    def _key(self, x, y, z):
        c = self.cell
        return (int(math.floor(x / c)), int(math.floor(y / c)), int(math.floor(z / c)))

    def nearest_dist(self, px, py, pz, search=2):
        """Find nearest distance, searching ±search cells in each direction."""
        cx, cy, cz = self._key(px, py, pz)
        best = float('inf')
        for dx in range(-search, search + 1):
            for dy in range(-search, search + 1):
                for dz in range(-search, search + 1):
                    for idx in self.grid.get((cx + dx, cy + dy, cz + dz), ()):
                        x, y, z = self.pts[idx]
                        d = (x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2
                        if d < best:
                            best = d
        return math.sqrt(best)


# ── Bounds helpers ────────────────────────────────────────────────────────────

def bounds(pts):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    zs = [p[2] for p in pts]
    return min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)


# ── Main validation ───────────────────────────────────────────────────────────

def validate(binary, cloud, mesh, expected_trees, min_coverage, no_run):
    stem = re.sub(r'\.ply$', '', cloud)
    seg_ply   = stem + '_segmented.ply'
    trees_mesh = stem + '_trees_mesh.ply'

    ok = True

    # ── 1. Run extraction ────────────────────────────────────────────────────
    if not no_run:
        print(f"Running: {binary} trees {cloud} {mesh}")
        result = subprocess.run(
            [binary, 'trees', cloud, mesh],
            capture_output=True, text=True
        )
        stdout = result.stdout + result.stderr
        print(stdout.rstrip())
        if result.returncode != 0:
            print(f"FAIL  extraction exited with code {result.returncode}")
            return False
    else:
        stdout = ''

    # ── 2. Tree count ────────────────────────────────────────────────────────
    if not no_run:
        m = re.search(r'(\d+) trees saved', stdout)
        if not m:
            print("FAIL  could not find 'N trees saved' in output")
            return False
        actual_trees = int(m.group(1))
        status = 'PASS' if actual_trees == expected_trees else 'FAIL'
        if status == 'FAIL':
            ok = False
        print(f"{status}  tree count: expected {expected_trees}, got {actual_trees}")

    # ── 3. Point retention ───────────────────────────────────────────────────
    if os.path.exists(seg_ply):
        n_input = read_ply_vertex_count(cloud)
        n_seg   = read_ply_vertex_count(seg_ply)
        retention = n_seg / n_input if n_input > 0 else 0.0
        status = 'PASS' if retention >= 0.95 else ('WARN' if retention >= 0.70 else 'FAIL')
        if status == 'FAIL':
            ok = False
        print(f"{status}  point retention: {n_seg}/{n_input} = {retention:.1%}")
    else:
        print(f"WARN  segmented output not found: {seg_ply}")

    # ── 4. Load mesh vertices for spatial checks ─────────────────────────────
    if not os.path.exists(trees_mesh):
        print(f"FAIL  tree mesh not found: {trees_mesh}")
        return False

    print(f"\nLoading mesh vertices from {os.path.basename(trees_mesh)} …")
    mesh_pts, n_mesh = read_ply_xyz(trees_mesh, sample_every=1)
    print(f"  {n_mesh} mesh vertices")

    # ── 5. Sample input points ───────────────────────────────────────────────
    # Every 200th point keeps runtime under a few seconds
    sample_rate = 200
    print(f"Sampling input cloud {os.path.basename(cloud)} (every {sample_rate}th point) …")
    sample_pts, n_cloud = read_ply_xyz(cloud, sample_every=sample_rate)
    print(f"  {len(sample_pts)} sampled points from {n_cloud} total")

    cloud_bnd  = bounds(sample_pts)
    mesh_bnd   = bounds(mesh_pts)

    # ── 6. Height completeness ───────────────────────────────────────────────
    cloud_zmax = cloud_bnd[5]
    mesh_zmax  = mesh_bnd[5]
    height_ratio = mesh_zmax / cloud_zmax if cloud_zmax > 0 else 0
    status = 'PASS' if height_ratio >= 0.85 else ('WARN' if height_ratio >= 0.70 else 'FAIL')
    if status == 'FAIL':
        ok = False
    print(f"\n{status}  height completeness: mesh_zmax={mesh_zmax:.2f}  cloud_zmax={cloud_zmax:.2f}  ratio={height_ratio:.2%}")

    # ── 7. XY extent ─────────────────────────────────────────────────────────
    cloud_xspan = cloud_bnd[1] - cloud_bnd[0]
    cloud_yspan = cloud_bnd[3] - cloud_bnd[2]
    mesh_xspan  = mesh_bnd[1] - mesh_bnd[0]
    mesh_yspan  = mesh_bnd[3] - mesh_bnd[2]
    x_ratio = mesh_xspan / cloud_xspan if cloud_xspan > 0 else 0
    y_ratio = mesh_yspan / cloud_yspan if cloud_yspan > 0 else 0
    xy_ratio = min(x_ratio, y_ratio)
    status = 'PASS' if xy_ratio >= 0.70 else ('WARN' if xy_ratio >= 0.50 else 'FAIL')
    if status == 'FAIL':
        ok = False
    print(f"{status}  XY coverage:  X={x_ratio:.0%}  Y={y_ratio:.0%}  (min={xy_ratio:.0%})")

    # ── 8. Point-to-mesh coverage ─────────────────────────────────────────────
    # For each sampled input point, find the distance to the nearest mesh vertex.
    # The mesh represents branch SKELETONS (cylinders), not the full crown volume,
    # so leaf/tip points can be several metres from the nearest cylinder surface.
    # Primary check: ≥ min_coverage of points within 2 m of some mesh vertex.
    # Points with no mesh vertex within the search range are reported separately.
    print(f"\nBuilding voxel grid on mesh …")
    grid = VoxelGrid(mesh_pts, cell_size=2.0)

    thresholds = [1.0, 2.0, 4.0]
    counts = [0] * len(thresholds)
    total = len(sample_pts)
    finite_dists = []
    n_unreachable = 0

    for px, py, pz in sample_pts:
        d = grid.nearest_dist(px, py, pz)
        if math.isinf(d):
            n_unreachable += 1
        else:
            finite_dists.append(d)
            for k, t in enumerate(thresholds):
                if d <= t:
                    counts[k] += 1

    if finite_dists:
        mean_d = sum(finite_dists) / len(finite_dists)
        finite_dists.sort()
        p50 = finite_dists[len(finite_dists) // 2]
        p95 = finite_dists[int(len(finite_dists) * 0.95)]
    else:
        mean_d = p50 = p95 = float('inf')

    print(f"  mean dist to nearest mesh vertex: {mean_d:.3f} m")
    print(f"  median: {p50:.3f} m   p95: {p95:.3f} m")
    if n_unreachable:
        print(f"  points with no mesh vertex within search range: {n_unreachable}/{total}")
    print()
    for k, t in enumerate(thresholds):
        cov = counts[k] / total
        is_key = (t == 2.0)
        status = 'PASS' if cov >= min_coverage else ('WARN' if cov >= min_coverage * 0.8 else 'FAIL')
        if is_key and status == 'FAIL':
            ok = False
        print(f"  {status}  coverage within {t:.1f} m: {counts[k]}/{total} = {cov:.1%}{'  ← primary check' if is_key else ''}")

    # ── 9. Phantom branch check ───────────────────────────────────────────────
    # Mesh vertices that are far from all input points indicate phantom branches.
    print(f"\nChecking for phantom mesh vertices …")
    cloud_grid = VoxelGrid(sample_pts, cell_size=2.0)
    far_count = sum(1 for x, y, z in mesh_pts if cloud_grid.nearest_dist(x, y, z) > 4.0)
    phantom_frac = far_count / len(mesh_pts) if mesh_pts else 1.0
    status = 'PASS' if phantom_frac < 0.10 else ('WARN' if phantom_frac < 0.25 else 'FAIL')
    if status == 'FAIL':
        ok = False
    print(f"  {status}  phantom vertices (>4 m from cloud): {far_count}/{len(mesh_pts)} = {phantom_frac:.1%}")

    # ── Summary ──────────────────────────────────────────────────────────────
    print()
    if ok:
        print("RESULT: PASS")
    else:
        print("RESULT: FAIL")
    return ok


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--binary', default='rayextract',
                   help='path to rayextract binary')
    p.add_argument('--cloud',  required=True,
                   help='input point cloud (.ply)')
    p.add_argument('--mesh',   required=True,
                   help='ground mesh (.ply)')
    p.add_argument('--expected-trees', type=int, default=1,
                   help='expected number of trees (default: 1)')
    p.add_argument('--min-coverage', type=float, default=0.70,
                   help='minimum fraction of points within 2 m of mesh (default: 0.70)')
    p.add_argument('--no-run', action='store_true',
                   help='skip extraction, validate existing output files')
    args = p.parse_args()

    ok = validate(
        binary=args.binary,
        cloud=os.path.abspath(args.cloud),
        mesh=os.path.abspath(args.mesh),
        expected_trees=args.expected_trees,
        min_coverage=args.min_coverage,
        no_run=args.no_run,
    )
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
