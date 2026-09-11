from pathlib import Path
import csv
import math

import numpy as np


FINAL_TIME = 0.01

CASES = {
    3000: Path(
        "results/convergence/time/exact/nt3000.csv"
    ),
    5999: Path(
        "results/convergence/time/exact/nt5999.csv"
    ),
    11997: Path(
        "results/convergence/time/exact/nt11997.csv"
    ),
}

FIELDS = [
    "T",
    "Y_CH4",
    "Y_O2",
    "Y_H2O",
    "Y_CO2",
]


def load_csv(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Missing convergence result: {path}"
        )

    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))

    return {
        field: np.asarray(
            [
                float(row[field])
                for row in rows
            ],
            dtype=float,
        )
        for field in FIELDS
    }


def l2_difference(a, b):
    return np.linalg.norm(
        a - b
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


def main():
    data = {
        nt: load_csv(path)
        for nt, path in CASES.items()
    }

    coarse_nt = 3000
    medium_nt = 5999
    fine_nt = 11997

    print()
    print("Temporal convergence study")
    print("==========================")
    print()

    for nt in [
        coarse_nt,
        medium_nt,
        fine_nt,
    ]:
        dt = (
            FINAL_TIME
            / (nt - 1)
        )

        tmax = np.max(
            data[nt]["T"]
        )

        print(
            f"nt={nt:5d}  "
            f"dt={dt:.10e} s  "
            f"Tmax={tmax:.9f} K"
        )

    print()
    print("Full-field convergence")
    print("======================")

    for field in FIELDS:
        coarse = data[coarse_nt][field]
        medium = data[medium_nt][field]
        fine = data[fine_nt][field]

        e_cm = l2_difference(
            coarse,
            medium,
        )

        e_mf = l2_difference(
            medium,
            fine,
        )

        rel_cm = relative_l2(
            coarse,
            medium,
        )

        rel_mf = relative_l2(
            medium,
            fine,
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
            f"  coarse -> medium relL2 : "
            f"{rel_cm:.8e}"
        )
        print(
            f"  medium -> fine   relL2 : "
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

    print()
    print("Tmax convergence")
    print("================")

    q1 = np.max(
        data[coarse_nt]["T"]
    )

    q2 = np.max(
        data[medium_nt]["T"]
    )

    q3 = np.max(
        data[fine_nt]["T"]
    )

    d12 = abs(
        q1 - q2
    )

    d23 = abs(
        q2 - q3
    )

    order = observed_order(
        d12,
        d23,
    )

    print(
        f"|T_coarse - T_medium| = "
        f"{d12:.9e} K"
    )

    print(
        f"|T_medium - T_fine|   = "
        f"{d23:.9e} K"
    )

    print(
        f"Observed Tmax order    = "
        f"{order:.6f}"
    )

    if (
        math.isfinite(order)
        and abs(2.0**order - 1.0) > 1.0e-14
    ):
        extrapolated = (
            q3
            + (q3 - q2)
            / (2.0**order - 1.0)
        )

        fine_error_estimate = abs(
            extrapolated - q3
        )

        print(
            f"Richardson Tmax(dt->0) = "
            f"{extrapolated:.9f} K"
        )

        print(
            f"Estimated fine-grid error = "
            f"{fine_error_estimate:.9e} K"
        )


if __name__ == "__main__":
    main()
