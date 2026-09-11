from pathlib import Path
import csv

import numpy as np


FILE = Path(
    "results/state_2999_reference_open_outlet.csv"
)

NX = 50
NY = 50

LX = 0.002
LY = 0.002

DX = LX / (NX - 1)
DY = LY / (NY - 1)


def load_field(rows, name):
    return np.asarray(
        [float(row[name]) for row in rows],
        dtype=float,
    ).reshape(NY, NX)


with FILE.open(newline="") as f:
    rows = list(csv.DictReader(f))

if len(rows) != NX * NY:
    raise RuntimeError(
        f"Expected {NX*NY} rows, got {len(rows)}"
    )

u = load_field(rows, "u")
v = load_field(rows, "v")


# ============================================================
# 1. Physical boundary-flux integral
# ============================================================

left_flux = -np.trapezoid(
    u[:, 0],
    dx=DY,
)

right_flux = np.trapezoid(
    u[:, -1],
    dx=DY,
)

bottom_flux = -np.trapezoid(
    v[0, :],
    dx=DX,
)

top_flux = np.trapezoid(
    v[-1, :],
    dx=DX,
)

physical_net = (
    left_flux
    + right_flux
    + bottom_flux
    + top_flux
)

physical_scale = (
    abs(left_flux)
    + abs(right_flux)
    + abs(bottom_flux)
    + abs(top_flux)
)


# ============================================================
# 2. Divergence used by the C++ projection
#
# D-u =
#
#   [u(i,j)-u(i-1,j)] / dx
# + [v(i,j)-v(i,j-1)] / dy
#
# on interior nodes.
# ============================================================

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

div_rms = np.sqrt(
    np.mean(divergence**2)
)

div_max = np.max(
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

normalized_rms = (
    div_rms
    / (
        velocity_scale
        / min(DX, DY)
    )
)


# ============================================================
# 3. Discrete divergence theorem for D-
#
# Summing D-u over all interior nodes telescopes exactly.
# ============================================================

operator_flux_x = (
    DY
    * np.sum(
        u[1:-1, -2]
        - u[1:-1, 0]
    )
)

operator_flux_y = (
    DX
    * np.sum(
        v[-2, 1:-1]
        - v[0, 1:-1]
    )
)

operator_net = (
    operator_flux_x
    + operator_flux_y
)

divergence_integral = (
    np.sum(divergence)
    * DX
    * DY
)


# ============================================================
# 4. Difference between numerical and physical boundaries
# ============================================================

right_velocity_gap = np.max(
    np.abs(
        u[:, -1]
        - u[:, -2]
    )
)

top_velocity_gap = np.max(
    np.abs(
        v[-1, :]
        - v[-2, :]
    )
)

bottom_velocity_gap = np.max(
    np.abs(
        v[0, :]
        - v[1, :]
    )
)

top_flux_gap = np.trapezoid(
    v[-1, :]
    - v[-2, :],
    dx=DX,
)

bottom_flux_gap = np.trapezoid(
    v[0, :]
    - v[1, :],
    dx=DX,
)


print()
print("Continuity / boundary diagnostic")
print("================================")

print()
print("1. Physical boundary flux")
print("-------------------------")

print(
    f"left        : {left_flux:+.12e}"
)

print(
    f"right       : {right_flux:+.12e}"
)

print(
    f"bottom      : {bottom_flux:+.12e}"
)

print(
    f"top         : {top_flux:+.12e}"
)

print(
    f"net         : {physical_net:+.12e}"
)

print(
    "relative net: "
    f"{abs(physical_net) / physical_scale:.12e}"
)


print()
print("2. Projection divergence")
print("------------------------")

print(
    f"RMS D-u       : {div_rms:.12e} s^-1"
)

print(
    f"max |D-u|     : {div_max:.12e} s^-1"
)

print(
    f"normalized RMS: {normalized_rms:.12e}"
)


print()
print("3. D- discrete flux balance")
print("---------------------------")

print(
    f"x contribution : {operator_flux_x:+.12e}"
)

print(
    f"y contribution : {operator_flux_y:+.12e}"
)

print(
    f"operator net   : {operator_net:+.12e}"
)

print(
    f"sum(div)*dA   : {divergence_integral:+.12e}"
)

print(
    "telescoping mismatch: "
    f"{abs(operator_net-divergence_integral):.12e}"
)


print()
print("4. Boundary/operator mismatch")
print("-----------------------------")

print(
    "max |u(right)-u(right-1)| : "
    f"{right_velocity_gap:.12e}"
)

print(
    "max |v(top)-v(top-1)|     : "
    f"{top_velocity_gap:.12e}"
)

print(
    "max |v(bottom)-v(bottom+1)|: "
    f"{bottom_velocity_gap:.12e}"
)

print(
    f"top flux gap    : {top_flux_gap:+.12e}"
)

print(
    f"bottom flux gap : {bottom_flux_gap:+.12e}"
)
