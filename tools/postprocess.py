#!/usr/bin/env python3
"""Post-processing for mdsim XYZ trajectories.

Computes the radial distribution function g(r) and the mean-square
displacement (MSD) from a multi-frame XYZ file written by `mdsim -traj`.

Usage:
    python postprocess.py run.xyz [--box L] [--rdf] [--msd] [--bins N]
                                  [--rmax R] [--dr 0.05] [--first N]
"""
import argparse
import math
import sys


def read_xyz(path):
    """Yield (comment, [(symbol, x, y, z), ...]) frames from an XYZ file."""
    frames = []
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue
            n = int(line)
            comment = f.readline().strip()
            atoms = []
            for _ in range(n):
                parts = f.readline().split()
                atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
            frames.append((comment, atoms))
    return frames


def compute_rdf(frames, box, rmax, bins):
    if len(frames) < 2:
        return None, None
    ref = frames[0][1]
    n = len(ref)
    dr = rmax / bins
    hist = [0] * bins
    pairs = 0
    for _, atoms in frames:
        for i in range(n):
            for j in range(i + 1, n):
                dx = atoms[i][1] - atoms[j][1]
                dy = atoms[i][2] - atoms[j][2]
                dz = atoms[i][3] - atoms[j][3]
                dx -= box * round(dx / box)
                dy -= box * round(dy / box)
                dz -= box * round(dz / box)
                r = math.sqrt(dx * dx + dy * dy + dz * dz)
                if r < rmax:
                    hist[int(r / dr)] += 1
                pairs += 1
    rho = n / (box ** 3)
    g = []
    for b in range(bins):
        r = (b + 0.5) * dr
        shell = 4.0 / 3.0 * math.pi * ((r + dr) ** 3 - r ** 3)
        g.append(hist[b] / (pairs * rho * shell))
    return [(b * dr, gv) for b, gv in enumerate(g)], rho


def compute_msd(frames, box):
    if len(frames) < 2:
        return None
    ref = frames[0][1]
    n = len(ref)
    msd = []
    for _, atoms in frames:
        acc = 0.0
        for i in range(n):
            dx = atoms[i][1] - ref[i][1]
            dy = atoms[i][2] - ref[i][2]
            dz = atoms[i][3] - ref[i][3]
            dx -= box * round(dx / box)
            dy -= box * round(dy / box)
            dz -= box * round(dz / box)
            acc += dx * dx + dy * dy + dz * dz
        msd.append(acc / n)
    return msd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--box", type=float, default=None, help="box side length (A)")
    ap.add_argument("--rdf", action="store_true")
    ap.add_argument("--msd", action="store_true")
    ap.add_argument("--bins", type=int, default=200)
    ap.add_argument("--rmax", type=float, default=5.0)
    args = ap.parse_args()

    frames = read_xyz(args.path)
    if not frames:
        print("no frames found", file=sys.stderr)
        return 1

    # Infer box from the frame comment if not given.
    box = args.box
    if box is None:
        for comment, _ in frames:
            toks = comment.split()
            for k, tok in enumerate(toks):
                if tok == "box" and k + 1 < len(toks):
                    try:
                        box = float(toks[k + 1])
                    except ValueError:
                        pass
                elif tok.startswith("box="):
                    try:
                        box = float(tok.split("=")[1])
                    except ValueError:
                        pass
        if box is None:
            print("could not infer box; pass --box", file=sys.stderr)
            return 1

    print(f"frames={len(frames)} atoms={len(frames[0][1])} box={box:.4f}")

    if args.rdf:
        g, rho = compute_rdf(frames, box, args.rmax, args.bins)
        print(f"# density = {rho:.6f} atoms/A^3")
        print("# r(A)  g(r)")
        for r, gv in g:
            print(f"{r:.4f} {gv:.6f}")

    if args.msd:
        msd = compute_msd(frames, box)
        print("# frame  MSD(A^2)")
        for k, m in enumerate(msd):
            print(f"{k} {m:.6f}")


if __name__ == "__main__":
    sys.exit(main())
