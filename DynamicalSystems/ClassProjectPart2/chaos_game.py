"""Reproduce Barnsley and Vince, Example 1; export data for PGFPlots."""
from pathlib import Path
import numpy as np
import argparse

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-only", action="store_true", help="Export CSV without Matplotlib")
    args = parser.parse_args()
    for folder in ("data", "plots"):
        (ROOT / folder).mkdir(exist_ok=True)
    rng = np.random.default_rng(535)
    steps = rng.integers(0, 2, size=4096)
    phases = (np.cumsum(steps) * ((np.sqrt(5) - 1) / 2)) % 1
    points = np.column_stack((np.cos(2 * np.pi * phases), np.sin(2 * np.pi * phases)))
    sizes = 2 ** np.arange(3, 13)
    distances = []
    for n in sizes:
        ordered = np.sort(np.unique(phases[:n]))
        gap = np.diff(np.r_[ordered, ordered[0] + 1]).max()
        distances.append(2 * np.sin(np.pi * gap / 2))
    assert np.all(np.diff(distances) <= 1e-14) and distances[-1] < 0.01
    np.savetxt(ROOT / "data/circle.csv", points[:512], delimiter=",", header="x,y", comments="")
    np.savetxt(ROOT / "data/coverage.csv", np.c_[sizes, distances], delimiter=",",
               header="n,distance", comments="")
    print(f"Coverage errors at N=8 and N=4096: {distances[0]:.6g}, {distances[-1]:.6g}")
    if args.data_only:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6), layout="constrained")
    axes[0].scatter(*points[:512].T, s=3)
    axes[0].set(aspect="equal", xlabel="x", ylabel="y", title="First 512 iterates")
    axes[1].loglog(sizes, distances, "o-")
    axes[1].set(xlabel="Number of iterates", ylabel="Hausdorff distance", title="Circle coverage")
    fig.savefig(ROOT / "plots/chaos_game.png", dpi=180)


if __name__ == "__main__":
    main()
