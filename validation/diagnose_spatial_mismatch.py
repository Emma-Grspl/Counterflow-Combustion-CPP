from pathlib import Path
import csv

import numpy as np


NX = 50

PYTHON_DATA = Path(
    "/Users/emma.grospellier/Projects/"
    "Final-project-M2-Combustion-of-methane-in-a-countercurrent-configuration/"
    "data"
)

CPP_FILE = Path(
    "results/state_2999_reference.csv"
)

FILES = {
    "Y_N2": "YN2.npy",
    "Y_CH4": "YCH4.npy",
    "Y_O2": "YO2.npy",
    "Y_H2O": "YH2O.npy",
    "Y_CO2": "YCO2.npy",
}

REACTIVE_FIELDS = [
    "Y_CH4",
    "Y_O2",
    "Y_H2O",
    "Y_CO2",
]


def load_python():
    result = {}

    for field, filename in FILES.items():
        array = np.load(
            PYTHON_DATA / filename
        )

        if array.shape != (NX * NX, 4):
            raise ValueError(
                f"{filename}: unexpected shape "
                f"{array.shape}"
            )

        result[field] = (
            array[:, -1]
            .reshape(NX, NX)
        )

    return result


def load_cpp():
    with CPP_FILE.open(newline="") as f:
        rows = list(csv.DictReader(f))

    if len(rows) != NX * NX:
        raise ValueError(
            f"Expected {NX*NX} C++ points, "
            f"found {len(rows)}"
        )

    result = {}

    for field in [
        "Y_N2",
        "Y_CH4",
        "Y_O2",
        "Y_H2O",
        "Y_CO2",
    ]:
        result[field] = np.asarray(
            [float(row[field]) for row in rows]
        ).reshape(NX, NX)

    return result


def rel_l2(reference, candidate):
    return (
        np.linalg.norm(candidate - reference)
        / np.linalg.norm(reference)
    )


def correlation(reference, candidate):
    return np.corrcoef(
        reference.ravel(),
        candidate.ravel(),
    )[0, 1]


def transforms(field):
    return {
        "identity": field,
        "flip_x": np.fliplr(field),
        "flip_y": np.flipud(field),
        "flip_xy": np.flipud(
            np.fliplr(field)
        ),
        "transpose": field.T,
        "transpose_flip_x": np.fliplr(
            field.T
        ),
        "transpose_flip_y": np.flipud(
            field.T
        ),
        "transpose_flip_xy": np.flipud(
            np.fliplr(field.T)
        ),
    }


def cropped(field, width):
    if width == 0:
        return field

    return field[
        width:-width,
        width:-width,
    ]


def main():
    python = load_python()
    cpp = load_cpp()

    print()
    print("Spatial mismatch diagnosis")
    print("==========================")

    print()
    print("1. Orientation check")
    print("--------------------")

    for field in REACTIVE_FIELDS:
        candidates = []

        for name, transformed in transforms(
            python[field]
        ).items():
            corr = correlation(
                transformed,
                cpp[field],
            )

            candidates.append(
                (corr, name)
            )

        candidates.sort(reverse=True)

        best_corr, best_name = (
            candidates[0]
        )

        identity_corr = correlation(
            python[field],
            cpp[field],
        )

        print(
            f"{field:7s} "
            f"identity={identity_corr:.6f}  "
            f"best={best_name:18s} "
            f"{best_corr:.6f}"
        )

    print()
    print("2. Boundary sensitivity")
    print("-----------------------")

    for field in REACTIVE_FIELDS:
        print()
        print(field)

        for width in [0, 1, 2, 3, 5]:
            py = cropped(
                python[field],
                width,
            )

            cc = cropped(
                cpp[field],
                width,
            )

            print(
                f"  remove {width:2d} cells: "
                f"relL2={rel_l2(py, cc):.6e}"
            )

    print()
    print("3. Extrema")
    print("----------")

    for field in REACTIVE_FIELDS:
        py_min = python[field].min()
        py_max = python[field].max()

        cc_min = cpp[field].min()
        cc_max = cpp[field].max()

        print(
            f"{field:7s} "
            f"Python=[{py_min:.6e}, "
            f"{py_max:.6e}]  "
            f"C++=[{cc_min:.6e}, "
            f"{cc_max:.6e}]"
        )

    print()
    print("4. Domain-integrated species")
    print("----------------------------")

    for field in REACTIVE_FIELDS:
        py_mean = python[field].mean()
        cc_mean = cpp[field].mean()

        relative_difference = (
            abs(cc_mean - py_mean)
            / abs(py_mean)
            if abs(py_mean) > 1.0e-30
            else np.nan
        )

        print(
            f"{field:7s} "
            f"Python mean={py_mean:.6e}  "
            f"C++ mean={cc_mean:.6e}  "
            f"rel.diff={relative_difference:.6%}"
        )

    print()
    print("5. Mixture closure")
    print("------------------")

    py_sum = (
        python["Y_N2"]
        + python["Y_CH4"]
        + python["Y_O2"]
        + python["Y_H2O"]
        + python["Y_CO2"]
    )

    cc_sum = (
        cpp["Y_N2"]
        + cpp["Y_CH4"]
        + cpp["Y_O2"]
        + cpp["Y_H2O"]
        + cpp["Y_CO2"]
    )

    print(
        "Python:"
    )
    print(
        f"  min sum(Y) = "
        f"{py_sum.min():.12e}"
    )
    print(
        f"  max sum(Y) = "
        f"{py_sum.max():.12e}"
    )
    print(
        f"  max error  = "
        f"{np.max(np.abs(py_sum - 1.0)):.12e}"
    )

    print(
        "C++:"
    )
    print(
        f"  min sum(Y) = "
        f"{cc_sum.min():.12e}"
    )
    print(
        f"  max sum(Y) = "
        f"{cc_sum.max():.12e}"
    )
    print(
        f"  max error  = "
        f"{np.max(np.abs(cc_sum - 1.0)):.12e}"
    )

    print()
    print("6. Product stoichiometry")
    print("------------------------")

    py_ratio = (
        python["Y_H2O"].max()
        / python["Y_CO2"].max()
    )

    cc_ratio = (
        cpp["Y_H2O"].max()
        / cpp["Y_CO2"].max()
    )

    theoretical_ratio = (
        2.0 * 18.01 / 44.00
    )

    print(
        f"Theoretical 2*MH2O/MCO2 : "
        f"{theoretical_ratio:.9f}"
    )

    print(
        f"Python                    : "
        f"{py_ratio:.9f}"
    )

    print(
        f"C++                       : "
        f"{cc_ratio:.9f}"
    )


if __name__ == "__main__":
    main()
