#include <cmath>
#include <iostream>

#include "counterflow/Grid.hpp"
#include "counterflow/NavierStokes.hpp"

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

    counterflow::Field2D u_old(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D v_old(
        grid.nx,
        grid.ny,
        0.0
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

    // Analytical test field:
    //
    // u(x,y) = x
    // v(x,y) = 0
    //
    // Therefore:
    // du/dx = 1
    // du/dy = 0
    // Laplacian(u) = 0

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            const std::size_t k = grid.index(i, j);

            u_old(i, j) = grid.x[k];
        }
    }

    const double dt = 0.1;
    const double nu = 0.01;

    counterflow::compute_intermediate_velocity(
        u_old,
        v_old,
        u_star,
        v_star,
        grid,
        dt,
        nu
    );

    // At the center:
    //
    // x = 0.5
    // u = 0.5
    // du/dx = 1
    //
    // u* = u - dt * u * du/dx
    //    = 0.5 - 0.1 * 0.5
    //    = 0.45

    if (!approximately_equal(
            u_star(2, 2),
            0.45,
            1.0e-12
        ))
    {
        std::cerr
            << "Unexpected u_star value: "
            << u_star(2, 2)
            << '\n';

        return 1;
    }

    if (!approximately_equal(
            v_star(2, 2),
            0.0,
            1.0e-12
        ))
    {
        std::cerr
            << "Unexpected v_star value.\n";

        return 1;
    }

    // Boundary values must be copied unchanged.

    if (!approximately_equal(
            u_star(0, 2),
            u_old(0, 2),
            1.0e-12
        ))
    {
        std::cerr
            << "Boundary velocity was not preserved.\n";

        return 1;
    }

    return 0;
}
