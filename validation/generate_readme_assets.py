from pathlib import Path
import csv
import math

import matplotlib.pyplot as plt
import numpy as np


# ============================================================
# Paths
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

REFERENCE_FILE = (
    ROOT
    / "results"
    / "state_2999_reference.csv"
)

TIME_FILES = {
    3000: (
        ROOT
        / "results/convergence/time/exact/nt3000.csv"
    ),
    5999: (
        ROOT
        / "results/convergence/time/exact/nt5999.csv"
    ),
    11997: (
        ROOT
        / "results/convergence/time/exact/nt11997.csv"
    ),
}

SPACE_FILES = {
    33: (
        ROOT
        / "results/convergence/space/n33.csv"
    ),
    65: (
        ROOT
        / "results/convergence/space/n65.csv"
    ),
    129: (
        ROOT
        / "results/convergence/space/n129.csv"
    ),
}

BENCHMARK_FILE = (
    ROOT
    / "benchmark"
    / "benchmark_summary.csv"
)

ASSET_DIR = ROOT / "assets"


# ============================================================
# Physical / reference values
# ============================================================

LX = 0.002
LY = 0.002

FINAL_TIME = 0.01

PYTHON_TMAX = 2323.01

PYTHON_PRODUCT_MAX = {
    "Y_H2O": 1.090003e-1,
    "Y_CO2": 1.331486e-1,
}

SPECIES = [
    "Y_CH4",
    "Y_O2",
    "Y_H2O",
    "Y_CO2",
]

CONVERGENCE_FIELDS = [
    "T",
    "Y_CH4",
    "Y_O2",
    "Y_H2O",
    "Y_CO2",
]


# ============================================================
# Plot configuration
# ============================================================

plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 220,
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})


# ============================================================
# Utilities
# ============================================================

def require_file(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Required file not found:\n{path}"
        )


def load_csv(path):
    require_file(path)

    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        raise ValueError(
            f"Empty CSV: {path}"
        )

    data = {}

    for key in rows[0].keys():
        data[key] = np.asarray(
            [float(row[key]) for row in rows],
            dtype=float,
        )

    return data


def load_grid(path):
    data = load_csv(path)

    order = np.lexsort(
        (
            data["x"],
            data["y"],
        )
    )

    data = {
        key: value[order]
        for key, value in data.items()
    }

    xs = np.unique(data["x"])
    ys = np.unique(data["y"])

    nx = len(xs)
    ny = len(ys)

    if nx * ny != len(data["x"]):
        raise ValueError(
            f"{path} does not contain a complete "
            "rectangular grid."
        )

    return data, xs, ys, nx, ny


def reshape(data, field, nx, ny):
    return data[field].reshape(
        ny,
        nx,
    )


def relative_l2(a, b):
    denominator = np.linalg.norm(b)

    if denominator <= 1.0e-30:
        return math.nan

    return (
        np.linalg.norm(a - b)
        / denominator
    )


def save_figure(fig, filename):
    path = ASSET_DIR / filename

    fig.savefig(
        path,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(fig)

    print(f"  wrote {path.relative_to(ROOT)}")


def add_panel_label(ax, label):
    ax.text(
        -0.10,
        1.05,
        label,
        transform=ax.transAxes,
        fontweight="bold",
        fontsize=13,
        va="top",
    )


# ============================================================
# Asset 1 — Final temperature
# ============================================================

def generate_temperature_asset():
    data, xs, ys, nx, ny = load_grid(
        REFERENCE_FILE
    )

    temperature = reshape(
        data,
        "T",
        nx,
        ny,
    )

    hotspot = int(
        np.argmax(data["T"])
    )

    hotspot_x = (
        data["x"][hotspot] * 1.0e3
    )

    hotspot_y = (
        data["y"][hotspot] * 1.0e3
    )

    tmax = data["T"][hotspot]

    extent = [
        xs[0] * 1.0e3,
        xs[-1] * 1.0e3,
        ys[0] * 1.0e3,
        ys[-1] * 1.0e3,
    ]

    fig, ax = plt.subplots(
        figsize=(7.2, 5.8)
    )

    image = ax.imshow(
        temperature,
        origin="lower",
        extent=extent,
        aspect="equal",
        cmap="inferno",
    )

    levels = np.linspace(
        float(np.min(temperature)),
        float(np.max(temperature)),
        9,
    )

    ax.contour(
        xs * 1.0e3,
        ys * 1.0e3,
        temperature,
        levels=levels,
        linewidths=0.45,
        alpha=0.55,
    )

    ax.scatter(
        [hotspot_x],
        [hotspot_y],
        marker="x",
        s=70,
        linewidths=2.0,
        label=(
            f"Hotspot: {tmax:.1f} K"
        ),
    )

    ax.set_title(
        "Final temperature field"
    )

    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")

    ax.legend(
        loc="upper right",
        frameon=True,
    )

    cbar = fig.colorbar(
        image,
        ax=ax,
    )

    cbar.set_label("Temperature [K]")

    fig.tight_layout()

    save_figure(
        fig,
        "temperature_field.png",
    )


# ============================================================
# Asset 2 — Species
# ============================================================

def generate_species_asset():
    data, xs, ys, nx, ny = load_grid(
        REFERENCE_FILE
    )

    labels = {
        "Y_CH4": r"$Y_{CH_4}$",
        "Y_O2": r"$Y_{O_2}$",
        "Y_H2O": r"$Y_{H_2O}$",
        "Y_CO2": r"$Y_{CO_2}$",
    }

    extent = [
        xs[0] * 1.0e3,
        xs[-1] * 1.0e3,
        ys[0] * 1.0e3,
        ys[-1] * 1.0e3,
    ]

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(10.5, 8.5),
        constrained_layout=True,
    )

    for ax, field in zip(
        axes.flat,
        SPECIES,
    ):
        values = reshape(
            data,
            field,
            nx,
            ny,
        )

        image = ax.imshow(
            values,
            origin="lower",
            extent=extent,
            aspect="equal",
            cmap="viridis",
        )

        ax.set_title(
            labels[field]
        )

        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")

        cbar = fig.colorbar(
            image,
            ax=ax,
            shrink=0.88,
        )

        cbar.set_label(
            "Mass fraction"
        )

    fig.suptitle(
        "Final reactive-species fields",
        fontsize=15,
    )

    save_figure(
        fig,
        "species_fields.png",
    )


# ============================================================
# Asset 3 — Flow
# ============================================================

def generate_flow_asset():
    data, xs, ys, nx, ny = load_grid(
        REFERENCE_FILE
    )

    u = reshape(
        data,
        "u",
        nx,
        ny,
    )

    v = reshape(
        data,
        "v",
        nx,
        ny,
    )

    pressure = reshape(
        data,
        "p",
        nx,
        ny,
    )

    speed = np.sqrt(
        u**2 + v**2
    )

    x_mm = xs * 1.0e3
    y_mm = ys * 1.0e3

    fig, ax = plt.subplots(
        figsize=(8.0, 5.8)
    )

    image = ax.imshow(
        speed,
        origin="lower",
        extent=[
            x_mm[0],
            x_mm[-1],
            y_mm[0],
            y_mm[-1],
        ],
        aspect="equal",
        cmap="viridis",
    )

    pressure_levels = np.linspace(
        float(np.min(pressure)),
        float(np.max(pressure)),
        10,
    )

    contours = ax.contour(
        x_mm,
        y_mm,
        pressure,
        levels=pressure_levels,
        linewidths=0.55,
        alpha=0.55,
    )

    ax.clabel(
        contours,
        inline=True,
        fontsize=7,
        fmt="%.2f",
    )

    ax.streamplot(
        x_mm,
        y_mm,
        u,
        v,
        density=1.05,
        linewidth=0.65,
        arrowsize=0.75,
    )

    ax.set_title(
        "Counterflow velocity field and pressure contours"
    )

    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")

    cbar = fig.colorbar(
        image,
        ax=ax,
    )

    cbar.set_label(
        r"Velocity magnitude [m s$^{-1}$]"
    )

    fig.tight_layout()

    save_figure(
        fig,
        "flow_field.png",
    )


# ============================================================
# Convergence helpers
# ============================================================

def get_flat_fields(path, fields):
    data = load_csv(path)

    return {
        field: data[field]
        for field in fields
    }


def temporal_errors():
    cases = {
        nt: get_flat_fields(
            path,
            CONVERGENCE_FIELDS,
        )
        for nt, path in TIME_FILES.items()
    }

    nts = [
        3000,
        5999,
        11997,
    ]

    dts = [
        FINAL_TIME / (nt - 1)
        for nt in nts
    ]

    result = {}

    for field in CONVERGENCE_FIELDS:
        result[field] = [
            relative_l2(
                cases[3000][field],
                cases[5999][field],
            ),
            relative_l2(
                cases[5999][field],
                cases[11997][field],
            ),
        ]

    return np.asarray(
        dts[:2],
        dtype=float,
    ), result


def spatial_errors():
    grids = {}

    for n, path in SPACE_FILES.items():
        data, _, _, nx, ny = load_grid(
            path
        )

        grids[n] = {
            field: reshape(
                data,
                field,
                nx,
                ny,
            )
            for field in CONVERGENCE_FIELDS
        }

    result = {}

    for field in CONVERGENCE_FIELDS:
        coarse = grids[33][field]
        medium = grids[65][field]
        fine = grids[129][field]

        medium_on_coarse = (
            medium[::2, ::2]
        )

        fine_on_medium = (
            fine[::2, ::2]
        )

        result[field] = [
            relative_l2(
                coarse,
                medium_on_coarse,
            ),
            relative_l2(
                medium,
                fine_on_medium,
            ),
        ]

    hs = np.asarray([
        LX / (33 - 1),
        LX / (65 - 1),
    ])

    return hs, result


# ============================================================
# Asset 4 — Convergence
# ============================================================

def generate_convergence_asset():
    dts, time_error = temporal_errors()
    hs, space_error = spatial_errors()

    labels = {
        "T": "T",
        "Y_CH4": "CH4",
        "Y_O2": "O2",
        "Y_H2O": "H2O",
        "Y_CO2": "CO2",
    }

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12.0, 5.0),
        constrained_layout=True,
    )

    # --------------------------------------------------------
    # Temporal
    # --------------------------------------------------------

    ax = axes[0]

    for field in CONVERGENCE_FIELDS:
        ax.loglog(
            dts,
            time_error[field],
            marker="o",
            linewidth=1.6,
            label=labels[field],
        )

    reference_time = (
        time_error["T"][0]
        * dts
        / dts[0]
    )

    ax.loglog(
        dts,
        reference_time,
        linestyle="--",
        linewidth=1.2,
        label=r"$O(\Delta t)$",
    )

    ax.set_title(
        "Temporal convergence"
    )

    ax.set_xlabel(
        r"Time step $\Delta t$ [s]"
    )

    ax.set_ylabel(
        "Relative L2 difference"
    )

    ax.grid(
        True,
        which="both",
        alpha=0.25,
    )

    ax.legend()

    ax.invert_xaxis()

    add_panel_label(
        ax,
        "(a)",
    )

    # --------------------------------------------------------
    # Spatial
    # --------------------------------------------------------

    ax = axes[1]

    for field in CONVERGENCE_FIELDS:
        ax.loglog(
            hs,
            space_error[field],
            marker="o",
            linewidth=1.6,
            label=labels[field],
        )

    reference_space = (
        space_error["T"][0]
        * hs
        / hs[0]
    )

    ax.loglog(
        hs,
        reference_space,
        linestyle="--",
        linewidth=1.2,
        label=r"$O(h)$",
    )

    ax.set_title(
        "Spatial convergence"
    )

    ax.set_xlabel(
        r"Grid spacing $h$ [m]"
    )

    ax.set_ylabel(
        "Relative L2 difference"
    )

    ax.grid(
        True,
        which="both",
        alpha=0.25,
    )

    ax.legend()

    ax.invert_xaxis()

    add_panel_label(
        ax,
        "(b)",
    )

    fig.suptitle(
        "Numerical convergence of the C++ solver",
        fontsize=15,
    )

    save_figure(
        fig,
        "convergence.png",
    )


# ============================================================
# Asset 5 — Python / C++ validation and precision
# ============================================================

def generate_validation_asset():
    data, _, _, nx, ny = load_grid(
        REFERENCE_FILE
    )

    cpp_tmax = float(
        np.max(data["T"])
    )

    cpp_h2o = float(
        np.max(data["Y_H2O"])
    )

    cpp_co2 = float(
        np.max(data["Y_CO2"])
    )

    relative_errors = np.asarray([
        abs(
            cpp_tmax
            - PYTHON_TMAX
        )
        / PYTHON_TMAX
        * 100.0,

        abs(
            cpp_h2o
            - PYTHON_PRODUCT_MAX["Y_H2O"]
        )
        / PYTHON_PRODUCT_MAX["Y_H2O"]
        * 100.0,

        abs(
            cpp_co2
            - PYTHON_PRODUCT_MAX["Y_CO2"]
        )
        / PYTHON_PRODUCT_MAX["Y_CO2"]
        * 100.0,
    ])

    mass_sum = (
        data["Y_N2"]
        + data["Y_CH4"]
        + data["Y_O2"]
        + data["Y_H2O"]
        + data["Y_CO2"]
    )

    mass_error = float(
        np.max(
            np.abs(
                mass_sum - 1.0
            )
        )
    )

    stoich_ratio = (
        2.0 * 18.01 / 44.00
    )

    stoich_error = float(
        np.max(
            np.abs(
                data["Y_H2O"]
                - stoich_ratio
                * data["Y_CO2"]
            )
        )
    )

    u = reshape(
        data,
        "u",
        nx,
        ny,
    )

    v = reshape(
        data,
        "v",
        nx,
        ny,
    )

    dx = LX / (nx - 1)
    dy = LY / (ny - 1)

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

    rms_divergence = float(
        np.sqrt(
            np.mean(
                divergence**2
            )
        )
    )

    velocity_scale = max(
        float(
            np.max(
                np.sqrt(
                    u**2 + v**2
                )
            )
        ),
        1.0e-30,
    )

    normalized_divergence = (
        rms_divergence
        / (
            velocity_scale
            / min(dx, dy)
        )
    )

    invariant_values = np.asarray([
        mass_error,
        stoich_error,
        normalized_divergence,
    ])

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(11.5, 4.8),
        constrained_layout=True,
    )

    # --------------------------------------------------------
    # Legacy benchmark error
    # --------------------------------------------------------

    ax = axes[0]

    names = [
        "Tmax",
        "H2O max",
        "CO2 max",
    ]

    bars = ax.bar(
        names,
        relative_errors,
    )

    ax.set_title(
        "Legacy Python vs C++"
    )

    ax.set_ylabel(
        "Relative difference [%]"
    )

    ax.set_ylim(
        0.0,
        max(relative_errors) * 1.5,
    )

    ax.grid(
        True,
        axis="y",
        alpha=0.25,
    )

    for bar, value in zip(
        bars,
        relative_errors,
    ):
        ax.text(
            bar.get_x()
            + bar.get_width() / 2.0,
            bar.get_height(),
            f"{value:.3f}%",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    add_panel_label(
        ax,
        "(a)",
    )

    # --------------------------------------------------------
    # Numerical invariants
    # --------------------------------------------------------

    ax = axes[1]

    invariant_names = [
        r"$|\Sigma Y-1|_{\max}$",
        "Stoich. residual",
        "Normalized div.",
    ]

    bars = ax.bar(
        invariant_names,
        invariant_values,
    )

    ax.set_yscale("log")

    ax.set_title(
        "C++ numerical precision"
    )

    ax.set_ylabel(
        "Residual"
    )

    ax.grid(
        True,
        axis="y",
        which="both",
        alpha=0.25,
    )

    for bar, value in zip(
        bars,
        invariant_values,
    ):
        ax.text(
            bar.get_x()
            + bar.get_width() / 2.0,
            value * 1.5,
            f"{value:.1e}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    add_panel_label(
        ax,
        "(b)",
    )

    fig.suptitle(
        "Reference agreement and numerical consistency",
        fontsize=15,
    )

    save_figure(
        fig,
        "python_cpp_validation.png",
    )


# ============================================================
# Asset 6 — Performance
# ============================================================

def load_benchmark_summary():
    require_file(
        BENCHMARK_FILE
    )

    result = {}

    with BENCHMARK_FILE.open(
        newline=""
    ) as f:
        reader = csv.DictReader(f)

        for row in reader:
            result[
                row["metric"]
            ] = float(
                row["value"]
            )

    return result


def generate_performance_asset():
    benchmark = load_benchmark_summary()

    python_pressure_us = (
        benchmark[
            "python_dense_pressure_median"
        ]
        * 1.0e6
    )

    cpp_pressure_us = (
        benchmark[
            "cpp_sparse_pressure_median"
        ]
        * 1.0e6
    )

    cpp_full = benchmark[
        "cpp_full_median"
    ]

    python_pressure_total = benchmark[
        "python_pressure_2999_estimate"
    ]

    cpp_pressure_total = benchmark[
        "cpp_pressure_2999_estimate"
    ]

    speedup = (
        python_pressure_us
        / cpp_pressure_us
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(11.8, 4.9),
        constrained_layout=True,
    )

    # --------------------------------------------------------
    # Per-solve latency
    # --------------------------------------------------------

    ax = axes[0]

    labels = [
        "Python dense",
        "C++ sparse",
    ]

    values = [
        python_pressure_us,
        cpp_pressure_us,
    ]

    bars = ax.barh(
        labels,
        values,
    )

    ax.set_xscale("log")

    ax.set_xlabel(
        "Median pressure-solve time [µs]"
    )

    ax.set_title(
        f"Pressure Poisson solve — {speedup:.0f}× speed-up"
    )

    ax.grid(
        True,
        axis="x",
        which="both",
        alpha=0.25,
    )

    for bar, value in zip(
        bars,
        values,
    ):
        ax.text(
            value * 1.12,
            bar.get_y()
            + bar.get_height() / 2.0,
            (
                f"{value / 1000.0:.1f} ms"
                if value >= 1000.0
                else f"{value:.1f} µs"
            ),
            va="center",
            fontsize=9,
        )

    add_panel_label(
        ax,
        "(a)",
    )

    # --------------------------------------------------------
    # End-to-end context
    # --------------------------------------------------------

    ax = axes[1]

    labels = [
        "Python pressure\n2999 solves",
        "C++ full\nsimulation",
        "C++ pressure\n2999 solves",
    ]

    values = [
        python_pressure_total,
        cpp_full,
        cpp_pressure_total,
    ]

    bars = ax.barh(
        labels,
        values,
    )

    ax.set_xscale("log")

    ax.set_xlabel(
        "Runtime [s]"
    )

    ax.set_title(
        "Runtime context"
    )

    ax.grid(
        True,
        axis="x",
        which="both",
        alpha=0.25,
    )

    for bar, value in zip(
        bars,
        values,
    ):
        ax.text(
            value * 1.12,
            bar.get_y()
            + bar.get_height() / 2.0,
            f"{value:.3g} s",
            va="center",
            fontsize=9,
        )

    add_panel_label(
        ax,
        "(b)",
    )

    fig.suptitle(
        "Performance impact of sparse pre-factorized pressure solves",
        fontsize=15,
    )

    save_figure(
        fig,
        "performance.png",
    )


# ============================================================
# Main
# ============================================================

def main():
    ASSET_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    required = [
        REFERENCE_FILE,
        *TIME_FILES.values(),
        *SPACE_FILES.values(),
        BENCHMARK_FILE,
    ]

    missing = [
        path
        for path in required
        if not path.exists()
    ]

    if missing:
        print(
            "Cannot generate README assets."
        )

        print(
            "Missing required files:"
        )

        for path in missing:
            print(
                f"  - {path.relative_to(ROOT)}"
            )

        raise SystemExit(1)

    print("Generating README assets")
    print("========================")

    generate_temperature_asset()
    generate_species_asset()
    generate_flow_asset()
    generate_convergence_asset()
    generate_validation_asset()
    generate_performance_asset()

    print()
    print(
        "README assets generated successfully."
    )


if __name__ == "__main__":
    main()
