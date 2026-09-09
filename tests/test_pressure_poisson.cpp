#include <cmath>
#include <iostream>

#include "counterflow/Grid.hpp"
#include "counterflow/PressurePoisson.hpp"

bool approximately_equal(double a, double b, double tolerance)
{
    return std::abs(a - b) < tolerance;
}

int main()
{
    const counterflow::Grid2D grid(
        1.0,
        1.0,
        5,
        5
    );

    // ========================================================
    // Test 1: pressure RHS
    // ========================================================

    counterflow::Field2D u_star(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D v_star(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D rhs(
        grid.nx,
        grid.ny,
        0.0
    );

    // u*(x,y) = x
    // v*(x,y) = y
    //
    // du*/dx = 1
    // dv*/dy = 1
    //
    // divergence = 2

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            const std::size_t k = grid.index(i, j);

            u_star(i, j) = grid.x[k];
            v_star(i, j) = grid.y[k];
        }
    }

    const double rho = 2.0;
    const double dt = 0.5;

    counterflow::compute_pressure_rhs(
        u_star,
        v_star,
        rhs,
        grid,
        rho,
        dt
    );

    // rho / dt = 4
    // divergence = 2
    //
    // rhs = 8

    if (!approximately_equal(
            rhs(2, 2),
            8.0,
            1.0e-12
        ))
    {
        std::cerr
            << "Unexpected pressure RHS value: "
            << rhs(2, 2)
            << '\n';

        return 1;
    }

    // ========================================================
    // Test 2: sparse Poisson solver
    // ========================================================

    counterflow::PressurePoissonSolver solver(grid);

    if (solver.system_size() != 25)
    {
        std::cerr
            << "Unexpected pressure system size.\n";

        return 1;
    }

    // A dense 25 x 25 matrix would contain 625 entries.
    // The sparse matrix should contain far fewer.

    if (solver.nonzero_count() >= 625)
    {
        std::cerr
            << "Pressure matrix is unexpectedly dense.\n";

        return 1;
    }

    counterflow::Field2D manufactured_rhs(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D exact_pressure(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D computed_pressure(
        grid.nx,
        grid.ny,
        0.0
    );

    // Manufactured discrete solution:
    //
    // p_i = dx^2 * [
    //     i(i - 1) - m(m - 1)
    // ]
    //
    // where m = nx - 1.
    //
    // This choice satisfies exactly:
    //
    // left:   p(1,j) - p(0,j) = 0
    // right:  p(nx-1,j) = 0
    // top/bottom: dp/dy = 0
    //
    // and its discrete second derivative in x is exactly 2.

    const double m =
        static_cast<double>(grid.nx - 1);

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            const double ii =
                static_cast<double>(i);

            exact_pressure(i, j) =
                grid.dx * grid.dx
                * (
                    ii * (ii - 1.0)
                    - m * (m - 1.0)
                );

            const bool interior =
                i > 0 &&
                i < grid.nx - 1 &&
                j > 0 &&
                j < grid.ny - 1;

            if (interior)
            {
                manufactured_rhs(i, j) = 2.0;
            }
        }
    }

    solver.solve(
        manufactured_rhs,
        computed_pressure
    );

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            if (!approximately_equal(
                    computed_pressure(i, j),
                    exact_pressure(i, j),
                    1.0e-10
                ))
            {
                std::cerr
                    << "Pressure mismatch at ("
                    << i
                    << ", "
                    << j
                    << "): expected "
                    << exact_pressure(i, j)
                    << ", got "
                    << computed_pressure(i, j)
                    << '\n';

                return 1;
            }
        }
    }

    return 0;
}
