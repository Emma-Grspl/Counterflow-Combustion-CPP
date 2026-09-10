from pathlib import Path
import csv
import math

import matplotlib.pyplot as plt
import numpy as np


NX = 50
NPOINTS = NX * NX

PYTHON_DATA = Path(
    "/Users/emma.grospellier/Projects/"
    "Final-project-M2-Combustion-of-methane-in-a-countercurrent-configuration/"
    "data"
)

CPP_FILE = Path(
    "results/state_2999_reference.csv"
)

OUTPUT_DIR = Path(
    "validation/results"
)

FIGURE_DIR = OUTPUT_DIR / "figures"

PYTHON_TMAX = 2323.01

PYTHON_FILES = {
    "Y_N2": "YN2.npy",
    "Y_CH4": "YCH4.npy",
    "Y_O2": "YO2.npy",
    "Y_H2O": "YH2O.npy",
    "Y_CO2": "YCO2.npy",
}

STRICT_FIELDS = [
    "Y_CH4",
    "Y_O2",
    "Y_H2O",
    "Y_CO2",
]


def load_python_field(path):
    array = np.load(path)

    if array.shape != (NPOINTS, 4):
        raise ValueError(
            f"{path.name}: expected shape "
            f"({NPOINTS}, 4), got {array.shape}"
        )

    # The fourth column is the final stored snapshot.
    return array[:, -1].astype(float)


def load_python_reference():
    fields = {}

    for field, filename in PYTHON_FILES.items():
        path = PYTHON_DATA / filename

        if not path.exists():
            raise FileNotFoundError(
                f"Missing Python reference: {path}"
            )

        fields[field] = load_python_field(path)

    return fields


def load_cpp_reference():
    if not CPP_FILE.exists():
        raise FileNotFoundError(
            f"Missing C++ reference: {CPP_FILE}"
        )

    with CPP_FILE.open(newline="") as f:
        rows = list(csv.DictReader(f))

    if len(rows) != NPOINTS:
        raise ValueError(
            f"Expected {NPOINTS} C++ points, "
            f"got {len(rows)}"
        )

    fields = {}

    for name in [
        "x",
        "y",
        "u",
        "v",
        "p",
        "Y_N2",
        "Y_CH4",
        "Y_O2",
        "Y_H2O",
        "Y_CO2",
        "T",
    ]:
        fields[name] = np.asarray(
            [float(row[name]) for row in rows],
            dtype=float,
        )

    return fields


def metrics(reference, candidate):
    error = candidate - reference

    ref_norm = np.linalg.norm(reference)

    rel_l2 = (
        np.linalg.norm(error) / ref_norm
        if ref_norm > 1.0e-30
        else math.nan
    )

    mae = np.mean(np.abs(error))

    rmse = np.sqrt(
        np.mean(error**2)
    )

    linf = np.max(
        np.abs(error)
    )

    if (
        np.std(reference) > 1.0e-30
        and np.std(candidate) > 1.0e-30
    ):
        correlation = np.corrcoef(
            reference,
            candidate,
        )[0, 1]
    else:
        correlation = math.nan

    return {
        "rel_l2": rel_l2,
        "mae": mae,
        "rmse": rmse,
        "linf": linf,
        "correlation": correlation,
    }


def save_comparison_figures(
    field,
    python_field,
    cpp_field,
):
    python_2d = python_field.reshape(
        NX,
        NX,
    )

    cpp_2d = cpp_field.reshape(
        NX,
        NX,
    )

    error_2d = np.abs(
        cpp_2d - python_2d
    )

    common_min = min(
        python_2d.min(),
        cpp_2d.min(),
    )

    common_max = max(
        python_2d.max(),
        cpp_2d.max(),
    )

    extent = [
        0.0,
        2.0,
        0.0,
        2.0,
    ]

    fig = plt.figure(
        figsize=(6, 5)
    )

    image = plt.imshow(
        python_2d,
        origin="lower",
        extent=extent,
        aspect="equal",
        vmin=common_min,
        vmax=common_max,
    )

    plt.xlabel("x [mm]")
    plt.ylabel("y [mm]")
    plt.title(
        f"Python reference — {field}"
    )

    plt.colorbar(image)
    plt.tight_layout()

    fig.savefig(
        FIGURE_DIR / f"{field}_python.png",
        dpi=180,
    )

    plt.close(fig)

    fig = plt.figure(
        figsize=(6, 5)
    )

    image = plt.imshow(
        cpp_2d,
        origin="lower",
        extent=extent,
        aspect="equal",
        vmin=common_min,
        vmax=common_max,
    )

    plt.xlabel("x [mm]")
    plt.ylabel("y [mm]")
    plt.title(
        f"C++ reference mode — {field}"
    )

    plt.colorbar(image)
    plt.tight_layout()

    fig.savefig(
        FIGURE_DIR / f"{field}_cpp.png",
        dpi=180,
    )

    plt.close(fig)

    fig = plt.figure(
        figsize=(6, 5)
    )

    image = plt.imshow(
        error_2d,
        origin="lower",
        extent=extent,
        aspect="equal",
    )

    plt.xlabel("x [mm]")
    plt.ylabel("y [mm]")
    plt.title(
        f"Absolute error — {field}"
    )

    plt.colorbar(image)
    plt.tight_layout()

    fig.savefig(
        FIGURE_DIR / f"{field}_error.png",
        dpi=180,
    )

    plt.close(fig)


def save_vertical_profile(
    field,
    python_field,
    cpp_field,
):
    python_2d = python_field.reshape(
        NX,
        NX,
    )

    cpp_2d = cpp_field.reshape(
        NX,
        NX,
    )

    # x ≈ 0.25 mm, inside the main inlet region.
    i = NX // 8

    y_mm = np.linspace(
        0.0,
        2.0,
        NX,
    )

    x_mm = (
        2.0 * i / (NX - 1)
    )

    fig = plt.figure(
        figsize=(7, 5)
    )

    plt.plot(
        y_mm,
        python_2d[:, i],
        label="Python",
    )

    plt.plot(
        y_mm,
        cpp_2d[:, i],
        "--",
        label="C++",
    )

    plt.xlabel("y [mm]")
    plt.ylabel(field)

    plt.title(
        f"{field} vertical profile "
        f"at x = {x_mm:.3f} mm"
    )

    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    fig.savefig(
        FIGURE_DIR
        / f"{field}_profile.png",
        dpi=180,
    )

    plt.close(fig)


def main():
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    FIGURE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    python = load_python_reference()
    cpp = load_cpp_reference()

    print()
    print("Python vs C++ spatial validation")
    print("================================")
    print()

    results = {}

    for field in STRICT_FIELDS:
        result = metrics(
            python[field],
            cpp[field],
        )

        results[field] = result

        print(
            f"{field:7s} "
            f"relL2={result['rel_l2']:.6e}  "
            f"RMSE={result['rmse']:.6e}  "
            f"Linf={result['linf']:.6e}  "
            f"corr={result['correlation']:.6f}"
        )

        save_comparison_figures(
            field,
            python[field],
            cpp[field],
        )

        save_vertical_profile(
            field,
            python[field],
            cpp[field],
        )

    print()
    print("Nitrogen diagnostic")
    print("===================")

    n2_result = metrics(
        python["Y_N2"],
        cpp["Y_N2"],
    )

    print(
        f"Y_N2    "
        f"relL2={n2_result['rel_l2']:.6e}  "
        f"corr={n2_result['correlation']:.6f}"
    )

    print()
    print("Temperature validation")
    print("======================")

    cpp_tmax = np.max(
        cpp["T"]
    )

    absolute_error = abs(
        cpp_tmax - PYTHON_TMAX
    )

    relative_error = (
        absolute_error / PYTHON_TMAX
    )

    hottest = int(
        np.argmax(cpp["T"])
    )

    print(
        f"Python Tmax : "
        f"{PYTHON_TMAX:.6f} K"
    )

    print(
        f"C++ Tmax    : "
        f"{cpp_tmax:.6f} K"
    )

    print(
        f"Abs error   : "
        f"{absolute_error:.6f} K"
    )

    print(
        f"Rel error   : "
        f"{relative_error:.6%}"
    )

    print(
        "C++ hotspot : "
        f"x={cpp['x'][hottest] * 1e3:.6f} mm, "
        f"y={cpp['y'][hottest] * 1e3:.6f} mm"
    )

    print()
    print("C++ mixture closure")
    print("===================")

    mass_sum = (
        cpp["Y_N2"]
        + cpp["Y_CH4"]
        + cpp["Y_O2"]
        + cpp["Y_H2O"]
        + cpp["Y_CO2"]
    )

    mass_error = np.max(
        np.abs(
            mass_sum - 1.0
        )
    )

    print(
        f"max |sum(Y)-1| = "
        f"{mass_error:.6e}"
    )

    metrics_path = (
        OUTPUT_DIR
        / "field_metrics.csv"
    )

    with metrics_path.open(
        "w",
        newline="",
    ) as f:
        writer = csv.writer(f)

        writer.writerow([
            "field",
            "relative_L2",
            "MAE",
            "RMSE",
            "Linf",
            "correlation",
        ])

        for field, result in results.items():
            writer.writerow([
                field,
                result["rel_l2"],
                result["mae"],
                result["rmse"],
                result["linf"],
                result["correlation"],
            ])

    print()
    print(
        f"Metrics: {metrics_path}"
    )

    print(
        f"Figures: {FIGURE_DIR}"
    )


if __name__ == "__main__":
    main()
