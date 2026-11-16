"""Sample the Sierpinski triangle with the chaos game and export plot data."""
import argparse
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
VERTICES = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-only", action="store_true", help="Export CSV without Matplotlib")
    args = parser.parse_args()
    for folder in ("data", "plots"):
        (ROOT / folder).mkdir(exist_ok=True)
    choices = np.random.default_rng(535).integers(0, 3, size=4096)
    points = np.empty((len(choices), 2))
    point = VERTICES[0].copy()
    for n, i in enumerate(choices):
        point = (point + VERTICES[i]) / 2
        points[n] = point

    np.savetxt(ROOT / "data/triangle.csv", points, delimiter=",",
               header="x,y", comments="", fmt="%.10f")
    if args.data_only:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6), layout="constrained")
    for ax, n in zip(axes, (256, len(points))):
        ax.scatter(*points[:n].T, s=2, color="tab:blue", linewidths=0)
        ax.set(aspect="equal", xlim=(-.03, 1.03), ylim=(-.03, .9),
               xlabel="x", ylabel="y", title=f"First {n:,} iterates")
    fig.savefig(ROOT / "plots/chaos_game.png", dpi=180)


if __name__ == "__main__":
    main()
