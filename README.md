# Counterflow Combustion C++ — Validated 2D Reactive-Flow Solver

A modern C++20 reimplementation and numerical reinterpretation of a Master's project on two-dimensional methane counterflow combustion. The project evolves the original Python prototype into a modular, validated C++ solver using sparse linear algebra, a compatible pressure-projection scheme, reactive advection-diffusion transport, adaptive chemistry subcycling, and systematic numerical verification. The original Python project is available in the GitHub repository: `Emma-Grspl/Final-project-M2-Combustion-of-methane-in-a-countercurrent-configuration`

---

## Highlights

- Modern C++20 implementation
- CMake-based build system
- Eigen sparse linear algebra
- Pre-factorized sparse pressure Poisson solver
- Compatible discrete pressure projection
- Explicit advection-diffusion transport
- Single-step methane oxidation chemistry
- Adaptive chemical subcycling
- Exact mixture mass-fraction closure
- Reference and automatic ignition modes
- Unit and numerical tests
- Temporal and spatial convergence studies
- Reproducible performance benchmarks

---

## Numerical Model

The solver considers a simplified two-dimensional methane counterflow configuration. The hydrodynamic component uses an incompressible fractional-step / projection method :

1. explicit velocity predictor,
2. pressure Poisson solve,
3. pressure-gradient correction,
4. physical velocity boundary conditions.

The discrete divergence and pressure-gradient operators are chosen to be compatible with the five-point Poisson stencil. Reactive scalars solve advection-diffusion-reaction equations for:

- methane, `CH4`,
- oxygen, `O2`,
- water, `H2O`,
- carbon dioxide, `CO2`,
- temperature.

Nitrogen is reconstructed from mixture closure:

$$Y_{N2} = 1 - Y_{CH4} - Y_{O2} - Y_{H2O} - Y_{CO2}$$

The chemistry uses a simplified single-step Arrhenius methane oxidation model:

$$CH_4 + 2 O_2 -> CO_2 + 2 H2O$$

Chemical evolution is subcycled inside each hydrodynamic time step.

The original Python project provided the physical case and reference solution, but the C++ version intentionally corrects several numerical issues. Notable changes include:

- sparse rather than dense pressure linear algebra,
- pressure-matrix factorization reused across time steps,
- discretely compatible pressure projection,
- explicit advection of reactive species,
- complete inlet species compositions,
- exact mixture mass-fraction closure,
- physically consistent open lateral outlet,
- modular separation between hydrodynamics, transport, chemistry, I/O, and validation.

Consequently, the C++ fields are not expected to match the legacy Python fields point-by-point. Validation therefore combines physical invariants, convergence studies, and selected global legacy-reference quantities.

---

## Validation

For the standard 50 x 50 reference configuration:

- Quantity	Result
- Maximum mixture-closure error	2.22e-16
- Maximum stoichiometric residual	5.07e-16
- Normalized RMS discrete divergence	7.43e-14
- Maximum advective CFL	9.61e-2
- 2-D diffusion number	6.00e-2
- Legacy Python benchmark
- Quantity	Python reference	C++	Relative difference
- Maximum temperature	2323.01 K	2320.95 K	0.0888 %
- Maximum H2O mass fraction	0.1090003	0.1088878	0.103 %
- Maximum CO2 mass fraction	0.1331486	0.1330111	0.103 %

The legacy comparison is used as a benchmark, not as a strict field-by-field regression test, because the C++ solver intentionally changes the discrete model.

---

## Convergence
### Temporal convergence

Keeping the spatial grid fixed and successively halving the time step gives approximately first-order convergence : 

- Temperature	0.998
- CH4 :	1.013
- O2 : 1.004
- H2O :	0.998
- CO2 :	0.998

The maximum-temperature diagnostic gives an observed order of approximately 0.988.

### Spatial convergence

Nested grids of 33 x 33, 65 x 65, and 129 x 129 give:

- u :	1.063
- v : 1.134
- p :	0.774
- T :	1.009
- CH4 :	1.135
- O2 : 0.870
- H2O	: 1.009
- CO2	: 1.009

The physical boundary-flux imbalance also decreases approximately as O(h) while the flux balance associated with the discrete projection operator remains at machine precision.

---

## Performance

Benchmarks were run in Release mode on the same local machine. End-to-end C++. For the standard 50 x 50, 2999-step simulation the median runtime is 13.55 s

### Pressure Poisson solve

| Metric                | Legacy Python |             C++ |
| --------------------- | ------------: | --------------: |
| Linear system size    |          2500 |            2500 |
| Representation        |         Dense |          Sparse |
| Matrix storage        |     ~47.7 MiB | 11,862 nonzeros |
| Median solve time     |      150.3 ms |        117.1 us |
| Estimated 2999 solves |       450.8 s |         0.351 s |

The sparse pre-factorized C++ pressure solve is approximately 1284x faster than the repeated dense NumPy solve measured in the legacy implementation.

This speed-up should not be interpreted as a pure C++ versus Python language comparison. It primarily reflects the numerical redesign from repeated dense factorizations to a reusable sparse factorization.

---

## Repository Structure
.
├── benchmark/      Performance benchmarks
├── include/
│   └── counterflow/
│                   Public C++ headers
├── src/            Solver implementation and executables
├── tests/          Unit and numerical tests
├── validation/     Validation and convergence utilities
├── .github/
│   └── workflows/  Continuous integration
├── CMakeLists.txt
├── LICENSE
├── README.md
└── requirements.txt

Generated simulation outputs and build artifacts are intentionally excluded from version control.

### Requirements

- C++ solver
- C++20-compatible compiler
- CMake >= 3.20
- Eigen >= 3.4

Tested with Apple Clang on macOS and automatically built/tested with GCC on Linux through GitHub Actions.

### Python analysis utilities

Python is only required for validation, plotting, and benchmarking utilities.

Create a virtual environment:

- python3 -m venv .venv
- source .venv/bin/activate
- python -m pip install --upgrade pip
- python -m pip install -r requirements.txt
- Build
- macOS with Homebrew

Install dependencies:

- brew install cmake eigen@3

Configure:

cmake \
  -S . \
  -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$(brew --prefix eigen@3)"

Build:

cmake --build build
Linux

With Eigen installed system-wide:

cmake \
  -S . \
  -B build \
  -DCMAKE_BUILD_TYPE=Release

cmake --build build
Tests

Run the C++ test suite with:

ctest \
  --test-dir build \
  --output-on-failure

The current suite contains nine test targets covering:

- configuration,
- grid and field storage,
- boundary conditions,
- Navier-Stokes projection,
- pressure Poisson solve,
- scalar transport,
- chemistry,
- reactive transport,
- complete simulation integration.
- Running the Solver

### Short run
./build/counterflow_combustion 100
Full reference run
./build/counterflow_combustion \
  2999 \
  results/state_2999_reference.csv \
  reference
Automatic steady-state activation
./build/counterflow_combustion \
  2999 \
  results/state_2999_auto.csv \
  auto

The reference mode reproduces the activation time used by the original Master's calculation. The auto mode activates thermal evolution after the C++ hydrodynamic steady-state criterion is satisfied.

---

## Numerical Validation

After generating the reference result: python validation/validate_solver.py

Temporal convergence: python validation/analyze_time_convergence.py

Spatial convergence: python validation/analyze_space_convergence.py

Boundary-flux convergence: python validation/analyze_boundary_flux_convergence.py

Additional diagnostic scripts are retained under validation/ to document differences between the legacy Python formulation and the C++ reinterpretation.

---

## Benchmarks

End-to-end C++ benchmark: python benchmark/benchmark_cpp.py

Pressure-solver microbenchmark: ./build/counterflow_pressure_benchmark

A summary of the measured reference results is stored in: benchmark/benchmark_summary.csv

---

## Scope and Limitations

This repository is a numerical-methods and scientific-computing project, not a production combustion CFD code.

The model intentionally uses:

- constant transport properties,
- simplified single-step methane chemistry,
- a two-dimensional geometry,
- explicit transport integration,
- a simplified incompressible formulation.

The project is primarily intended to demonstrate numerical modeling, modern C++ design, sparse linear algebra, verification, validation, and performance engineering.

---

## Original Project

This repository is derived from my Master's project: `Emma-Grspl/Final-project-M2-Combustion-of-methane-in-a-countercurrent-configuration`

The original implementation was written in Python and served as the physical and numerical starting point for this C++ reinterpretation.

## License

This project is released under the MIT License. See LICENSE for details.
