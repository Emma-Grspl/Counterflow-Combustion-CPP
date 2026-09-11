from pathlib import Path
import csv

import numpy as np


LX = 0.002
LY = 0.002

CASES = {
    33: Path(
        "results/convergence/space/n33.csv"
    ),
    65: Path(
        "results/convergence/space/n65.csv"
    ),
    129: Path(
        "results/convergence/space/n129.csv"
    ),
}


def load_velocity(path, n):
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))

    if len(rows) != n * n:
        raise RuntimeError(
            f"{path}: expected {n*n} rows, "
            f"found {len(rows)}"
        )

    u = np.asarray(
        [float(row["u"]) for row in rows],
        dtype=float,
    ).reshape(n, n)

    v = np.asarray(
        [float(row["v"]) for row in rows],
        dtype=float,
    ).reshape(n, n)

    return u, v


print()
print("Boundary-flux convergence")
print("=========================")
print()

previous_error = None

for n, path in CASES.items():
    u, v = load_velocity(
        path,
        n,
    )

    dx = LX / (n - 1)
    dy = LY / (n - 1)

    left = -np.trapezoid(
        u[:, 0],
        dx=dy,
    )

    right = np.trapezoid(
        u[:, -1],
        dx=dy,
    )

    bottom = -np.trapezoid(
        v[0, :],
        dx=dx,
    )

    top = np.trapezoid(
        v[-1, :],
        dx=dx,
    )

    net = (
        left
        + right
        + bottom
        + top
    )

    scale = (
        abs(left)
        + abs(right)
        + abs(bottom)
        + abs(top)
    )

    relative = (
        abs(net) / scale
    )

    # --------------------------------------------------------
    # D- divergence actually used by the projection
    # --------------------------------------------------------

    divergence = (
        (
            u[1:-1, 1:-1]
            - u[1:-1, :-2]
        )
        / dx
        +
        (
            v[1:-1, 1:-1]
            - v[:-2, 1:-1]
        )
        / dy
    )

    operator_balance = (
        np.sum(divergence)
        * dx
        * dy
    )

    if previous_error is None:
        observed_order = float("nan")
    else:
        observed_order = (
            np.log(
                previous_error
                / relative
            )
            / np.log(2.0)
        )

    print(
        f"N={n:3d}"
    )

    print(
        f"  physical net flux    : "
        f"{net:+.12e}"
    )

    print(
        f"  relative imbalance   : "
        f"{relative:.12e}"
    )

    print(
        f"  D- operator balance  : "
        f"{operator_balance:+.12e}"
    )

    if previous_error is not None:
        print(
            f"  imbalance order      : "
            f"{observed_order:.6f}"
        )

    print()

    previous_error = relative
