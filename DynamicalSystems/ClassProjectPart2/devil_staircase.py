"""Reproduce Bak's Figure 5c numerically and export data for PGFPlots."""
from pathlib import Path
import numpy as np
import argparse

ROOT = Path(__file__).resolve().parent


def winding(bias, burn, samples, coupling=1):
    phase = np.full_like(bias, 0.123)
    displacement = np.zeros_like(bias)
    for n in range(burn + samples):
        increment = bias + coupling / (2 * np.pi) * np.sin(2 * np.pi * phase)
        phase = (phase + increment) % 1
        if n >= burn:
            displacement += increment
    return displacement / samples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-only", action="store_true", help="Export CSV without Matplotlib")
    args = parser.parse_args()
    for folder in ("data", "plots"):
        (ROOT / folder).mkdir(exist_ok=True)
    bias = np.linspace(0, 1, 3001)
    w = winding(bias, 2048, 8192)
    w_long = winding(bias, 4096, 16384)
    assert np.allclose(winding(bias[:10], 4, 16, coupling=0), bias[:10])
    assert np.max(np.abs(w - w_long)) < 0.001
    for target in (0, 1/3, 1/2, 2/3, 1):
        locked = bias[np.abs(w_long - target) <= 2e-5]
        assert len(locked) > 1
        print(f"W={target:.6g}: grid interval [{locked[0]:.6f}, {locked[-1]:.6f}]")
    print(f"Maximum change under doubling: {np.max(np.abs(w - w_long)):.6g}")
    np.savetxt(ROOT / "data/staircase.csv", np.c_[bias, w_long], delimiter=",",
               header="bias,winding", comments="")
    if args.data_only:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 4), layout="constrained")
    ax.plot(bias, w_long, color="tab:blue")
    ax.plot(bias, bias, "--", color="0.5", label="Uncoupled")
    ax.set(xlabel="Bias", ylabel="Winding number", title="Critical sine circle map")
    ax.legend()
    fig.savefig(ROOT / "plots/devil_staircase.png", dpi=180)


if __name__ == "__main__":
    main()
