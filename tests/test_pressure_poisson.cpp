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

    // Analytical velocity field:
    //
    // u*(x,y) = x
    // v*(x,y) = y
    //
    // Therefore:
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
    // rhs = 4 * 2 = 8

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

    // Boundaries are not modified by compute_pressure_rhs.

    if (!approximately_equal(
            rhs(0, 0),
            0.0,
            1.0e-12
        ))
    {
        std::cerr
            << "Pressure RHS boundary was unexpectedly modified.\n";

        return 1;
    }

    return 0;
}
