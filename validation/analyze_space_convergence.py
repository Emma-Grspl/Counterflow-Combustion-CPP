from pathlib import Path
import csv
import math

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

FIELDS = [
    "u",
    "v",
    "p",
    "T",
    "Y_CH4",
    "Y_O2",
    "Y_H2O",
    "Y_CO2",
]


def load_csv(path, n):
    if not path.exists():
        raise FileNotFoundError(
            f"Missing spatial-convergence file: {path}"
        )

    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))

    expected = n * n

    if len(rows) != expected:
        raise ValueError(
            f"{path}: expected {expected} points, "
            f"found {len(rows)}"
        )

    result = {}

    for field in FIELDS:
        flat = np.asarray(
            [
                float(row[field])
                for row in rows
            ],
            dtype=float,
        )

        result[field] = flat.reshape(
            n,
            n,
        )

    return result


def rms_difference(a, b):
    return np.sqrt(
        np.mean(
            (a - b) ** 2
        )
    )


def relative_l2(a, b):
    denominator = np.linalg.norm(b)

    if denominator <= 1.0e-30:
        return math.nan

    return (
        np.linalg.norm(a - b)
        / denominator
    )


def observed_order(e_coarse_medium, e_medium_fine):
    if (
        e_coarse_medium <= 0.0
        or e_medium_fine <= 0.0
    ):
        return math.nan

    return (
        math.log(
            e_coarse_medium
            / e_medium_fine
        )
        / math.log(2.0)
    )


def restrict_fine_to_coarse(fine):
    """
    The grids are nested:

        33  -> 32 intervals
        65  -> 64 intervals
        129 -> 128 intervals

    Every second point of the fine grid therefore lies
    exactly on a point of the next coarser grid.
    """
    return fine[::2, ::2]


def main():
    data = {
        n: load_csv(path, n)
        for n, path in CASES.items()
    }

    print()
    print("Spatial convergence study")
    print("=========================")
    print()

    for n in [33, 65, 129]:
        dx = LX / (n - 1)
        dy = LY / (n - 1)

        print(
            f"N={n:3d} x {n:<3d}  "
            f"dx={dx:.10e} m  "
            f"dy={dy:.10e} m  "
            f"Tmax={np.max(data[n]['T']):.9f} K"
        )

    # --------------------------------------------------------
    # Nested-grid restriction
    # --------------------------------------------------------

    medium_on_coarse = {
        field: restrict_fine_to_coarse(
            data[65][field]
        )
        for field in FIELDS
    }

    fine_on_medium = {
        field: restrict_fine_to_coarse(
            data[129][field]
        )
        for field in FIELDS
    }

    print()
    print("Full-field spatial convergence")
    print("==============================")

    for field in FIELDS:
        coarse = data[33][field]
        medium = data[65][field]
        fine = data[129][field]

        medium_c = medium_on_coarse[field]
        fine_m = fine_on_medium[field]

        if coarse.shape != medium_c.shape:
            raise RuntimeError(
                f"{field}: coarse restriction mismatch "
                f"{coarse.shape} vs {medium_c.shape}"
            )

        if medium.shape != fine_m.shape:
            raise RuntimeError(
                f"{field}: medium restriction mismatch "
                f"{medium.shape} vs {fine_m.shape}"
            )

        e_cm = rms_difference(
            coarse,
            medium_c,
        )

        e_mf = rms_difference(
            medium,
            fine_m,
        )

        rel_cm = relative_l2(
            coarse,
            medium_c,
        )

        rel_mf = relative_l2(
            medium,
            fine_m,
        )

        order = observed_order(
            e_cm,
            e_mf,
        )

        reduction = (
            e_cm / e_mf
            if e_mf > 0.0
            else math.inf
        )

        print()
        print(field)

        print(
            f"  33 -> 65  RMS error    : "
            f"{e_cm:.8e}"
        )

        print(
            f"  65 -> 129 RMS error    : "
            f"{e_mf:.8e}"
        )

        print(
            f"  33 -> 65  relative L2  : "
            f"{rel_cm:.8e}"
        )

        print(
            f"  65 -> 129 relative L2  : "
            f"{rel_mf:.8e}"
        )

        print(
            f"  error reduction ratio  : "
            f"{reduction:.6f}"
        )

        print(
            f"  observed order p        : "
            f"{order:.6f}"
        )

    # --------------------------------------------------------
    # Integral quantities
    # --------------------------------------------------------

    print()
    print("Domain-mean convergence")
    print("=======================")

    for field in [
        "T",
        "Y_CH4",
        "Y_O2",
        "Y_H2O",
        "Y_CO2",
    ]:
        q33 = np.mean(
            data[33][field]
        )

        q65 = np.mean(
            data[65][field]
        )

        q129 = np.mean(
            data[129][field]
        )

        d1 = abs(
            q33 - q65
        )

        d2 = abs(
            q65 - q129
        )

        order = observed_order(
            d1,
            d2,
        )

        print(
            f"{field:7s} "
            f"N33={q33:.9e}  "
            f"N65={q65:.9e}  "
            f"N129={q129:.9e}  "
            f"p={order:.6f}"
        )

    # --------------------------------------------------------
    # Maximum temperature: diagnostic only
    # --------------------------------------------------------

    print()
    print("Tmax diagnostic")
    print("===============")

    t33 = np.max(
        data[33]["T"]
    )

    t65 = np.max(
        data[65]["T"]
    )

    t129 = np.max(
        data[129]["T"]
    )

    d1 = abs(
        t33 - t65
    )

    d2 = abs(
        t65 - t129
    )

    order = observed_order(
        d1,
        d2,
    )

    print(
        f"Tmax N=33  : {t33:.9f} K"
    )

    print(
        f"Tmax N=65  : {t65:.9f} K"
    )

    print(
        f"Tmax N=129 : {t129:.9f} K"
    )

    print(
        f"|T33-T65|  : {d1:.9e} K"
    )

    print(
        f"|T65-T129| : {d2:.9e} K"
    )

    print(
        f"Diagnostic order: {order:.6f}"
    )


if __name__ == "__main__":
    main()
