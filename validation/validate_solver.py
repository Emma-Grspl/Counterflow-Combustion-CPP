from pathlib import Path
import csv
import math
import sys

import numpy as np


CPP_FILE = Path(
    "results/state_2999_reference.csv"
)

# ------------------------------------------------------------
# Physical / numerical parameters
# ------------------------------------------------------------

NX = 50
NY = 50

LX = 0.002
LY = 0.002

FINAL_TIME = 0.01
NT = 3000

DIFFUSIVITY = 15.0e-6

DX = LX / (NX - 1)
DY = LY / (NY - 1)

DT = FINAL_TIME / (NT - 1)


# ------------------------------------------------------------
# Legacy Python benchmark values
# ------------------------------------------------------------

PYTHON_TMAX = 2323.01

PYTHON_PRODUCT_MAX = {
    "Y_H2O": 1.090003e-1,
    "Y_CO2": 1.331486e-1,
}


# ------------------------------------------------------------
# Acceptance tolerances
# ------------------------------------------------------------

MASS_CLOSURE_TOL = 1.0e-12

SPECIES_BOUND_TOL = 1.0e-12

NORMALIZED_DIVERGENCE_TOL = 1.0e-10

STOICHIOMETRY_TOL = 1.0e-10

TMAX_RELATIVE_TOL = 2.0e-2

PRODUCT_MAX_RELATIVE_TOL = 2.0e-2


SPECIES = [
    "Y_N2",
    "Y_CH4",
    "Y_O2",
    "Y_H2O",
    "Y_CO2",
]


def load_cpp_result(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Missing C++ result: {path}"
        )

    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))

    expected_size = NX * NY

    if len(rows) != expected_size:
        raise ValueError(
            f"Expected {expected_size} rows, "
            f"got {len(rows)}."
        )

    fields = {}

    for name in [
        "x",
        "y",
        "u",
        "v",
        "p",
        *SPECIES,
        "T",
    ]:
        fields[name] = np.asarray(
            [
                float(row[name])
                for row in rows
            ],
            dtype=float,
        )

    return fields


def relative_error(value, reference):
    return (
        abs(value - reference)
        / abs(reference)
    )


def print_check(
    name,
    value,
    tolerance,
    passed,
):
    status = (
        "PASS"
        if passed
        else "FAIL"
    )

    print(
        f"[{status}] "
        f"{name:<34s} "
        f"{value:.6e} "
        f"(tol={tolerance:.1e})"
    )


def main():
    data = load_cpp_result(
        CPP_FILE
    )

    failures = []

    print()
    print("Counterflow C++ solver validation")
    print("=================================")

    # ========================================================
    # 1. Finite values
    # ========================================================

    print()
    print("1. Finite-value checks")
    print("----------------------")

    finite_ok = True

    for field, values in data.items():
        if not np.all(
            np.isfinite(values)
        ):
            finite_ok = False
            failures.append(
                f"non-finite values in {field}"
            )

    print(
        "[PASS] All fields are finite"
        if finite_ok
        else "[FAIL] Non-finite values detected"
    )

    # ========================================================
    # 2. Species bounds
    # ========================================================

    print()
    print("2. Species bounds")
    print("-----------------")

    species_min = min(
        np.min(data[field])
        for field in SPECIES
    )

    species_max = max(
        np.max(data[field])
        for field in SPECIES
    )

    lower_violation = max(
        0.0,
        -species_min,
    )

    upper_violation = max(
        0.0,
        species_max - 1.0,
    )

    bound_violation = max(
        lower_violation,
        upper_violation,
    )

    bounds_ok = (
        bound_violation
        <= SPECIES_BOUND_TOL
    )

    print_check(
        "max species-bound violation",
        bound_violation,
        SPECIES_BOUND_TOL,
        bounds_ok,
    )

    if not bounds_ok:
        failures.append(
            "species outside [0,1]"
        )

    # ========================================================
    # 3. Mixture mass closure
    # ========================================================

    print()
    print("3. Mixture closure")
    print("------------------")

    mass_sum = sum(
        data[field]
        for field in SPECIES
    )

    mass_error = np.max(
        np.abs(
            mass_sum - 1.0
        )
    )

    closure_ok = (
        mass_error
        <= MASS_CLOSURE_TOL
    )

    print_check(
        "max |sum(Y)-1|",
        mass_error,
        MASS_CLOSURE_TOL,
        closure_ok,
    )

    if not closure_ok:
        failures.append(
            "mixture mass closure"
        )

    # ========================================================
    # 4. Reaction stoichiometry
    # ========================================================

    print()
    print("4. Product stoichiometry")
    print("------------------------")

    MH2O = 18.01
    MCO2 = 44.00

    expected_ratio = (
        2.0 * MH2O / MCO2
    )

    # With zero initial products and identical transport
    # operators, the product fields must preserve this
    # mass-stoichiometric relation.
    stoich_residual = np.max(
        np.abs(
            data["Y_H2O"]
            - expected_ratio
            * data["Y_CO2"]
        )
    )

    stoich_ok = (
        stoich_residual
        <= STOICHIOMETRY_TOL
    )

    print(
        f"Expected H2O/CO2 mass ratio: "
        f"{expected_ratio:.12f}"
    )

    print_check(
        "max stoichiometric residual",
        stoich_residual,
        STOICHIOMETRY_TOL,
        stoich_ok,
    )

    if not stoich_ok:
        failures.append(
            "product stoichiometry"
        )

    # ========================================================
    # 5. Discrete incompressibility
    # ========================================================

    print()
    print("5. Discrete incompressibility")
    print("-----------------------------")

    u = data["u"].reshape(
        NY,
        NX,
    )

    v = data["v"].reshape(
        NY,
        NX,
    )

    # This is the same backward-divergence operator D-
    # used by the pressure projection in the C++ solver.
    divergence = (
        (
            u[1:-1, 1:-1]
            - u[1:-1, :-2]
        )
        / DX
        +
        (
            v[1:-1, 1:-1]
            - v[:-2, 1:-1]
        )
        / DY
    )

    divergence_rms = np.sqrt(
        np.mean(
            divergence**2
        )
    )

    divergence_max = np.max(
        np.abs(divergence)
    )

    velocity_scale = max(
        np.max(
            np.sqrt(
                u**2 + v**2
            )
        ),
        1.0e-30,
    )

    divergence_scale = (
        velocity_scale
        / min(DX, DY)
    )

    normalized_divergence = (
        divergence_rms
        / divergence_scale
    )

    divergence_ok = (
        normalized_divergence
        <= NORMALIZED_DIVERGENCE_TOL
    )

    print(
        f"RMS D-u     : "
        f"{divergence_rms:.6e} s^-1"
    )

    print(
        f"max |D-u|   : "
        f"{divergence_max:.6e} s^-1"
    )

    print_check(
        "normalized RMS divergence",
        normalized_divergence,
        NORMALIZED_DIVERGENCE_TOL,
        divergence_ok,
    )

    if not divergence_ok:
        failures.append(
            "discrete incompressibility"
        )

    # ========================================================
    # 6. Explicit transport stability indicators
    # ========================================================

    print()
    print("6. Explicit transport indicators")
    print("--------------------------------")

    advective_cfl = np.max(
        np.abs(u) * DT / DX
        + np.abs(v) * DT / DY
    )

    diffusion_number = (
        DIFFUSIVITY
        * DT
        * (
            1.0 / DX**2
            + 1.0 / DY**2
        )
    )

    print(
        f"Max advective CFL : "
        f"{advective_cfl:.6e}"
    )

    print(
        f"2-D diffusion number : "
        f"{diffusion_number:.6e}"
    )

    # These are reported rather than used as strict
    # acceptance tests because chemistry is subcycled
    # independently.

    # ========================================================
    # 7. Legacy Python benchmark
    # ========================================================

    print()
    print("7. Legacy Python benchmark")
    print("--------------------------")

    cpp_tmax = np.max(
        data["T"]
    )

    tmax_relative_error = (
        relative_error(
            cpp_tmax,
            PYTHON_TMAX,
        )
    )

    tmax_ok = (
        tmax_relative_error
        <= TMAX_RELATIVE_TOL
    )

    print(
        f"Python Tmax : "
        f"{PYTHON_TMAX:.6f} K"
    )

    print(
        f"C++ Tmax    : "
        f"{cpp_tmax:.6f} K"
    )

    print_check(
        "Tmax relative error",
        tmax_relative_error,
        TMAX_RELATIVE_TOL,
        tmax_ok,
    )

    if not tmax_ok:
        failures.append(
            "legacy Tmax benchmark"
        )

    for field, reference_max in (
        PYTHON_PRODUCT_MAX.items()
    ):
        cpp_max = np.max(
            data[field]
        )

        error = relative_error(
            cpp_max,
            reference_max,
        )

        passed = (
            error
            <= PRODUCT_MAX_RELATIVE_TOL
        )

        print(
            f"{field} Python max : "
            f"{reference_max:.9e}"
        )

        print(
            f"{field} C++ max    : "
            f"{cpp_max:.9e}"
        )

        print_check(
            f"{field} max relative error",
            error,
            PRODUCT_MAX_RELATIVE_TOL,
            passed,
        )

        if not passed:
            failures.append(
                f"{field} legacy maximum"
            )

    # ========================================================
    # Final verdict
    # ========================================================

    print()
    print("=================================")

    if failures:
        print(
            "VALIDATION RESULT: FAIL"
        )

        print()
        print(
            "Failed checks:"
        )

        for failure in failures:
            print(
                f"  - {failure}"
            )

        return 1

    print(
        "VALIDATION RESULT: PASS"
    )

    print(
        "All strict physical and numerical "
        "checks passed."
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
