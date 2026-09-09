#include <cmath>
#include <iostream>

#include "counterflow/Grid.hpp"
#include "counterflow/Transport.hpp"

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
    // Test 1: nitrogen initialization
    // ========================================================

    counterflow::Field2D nitrogen(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::initialize_nitrogen_mass_fraction(
        nitrogen,
        grid
    );

    if (!approximately_equal(
            nitrogen(0, 0),
            0.79,
            1.0e-12
        ))
    {
        std::cerr
            << "Incorrect nitrogen inlet value.\n";

        return 1;
    }

    if (!approximately_equal(
            nitrogen(1, 0),
            1.0,
            1.0e-12
        ))
    {
        std::cerr
            << "Incorrect nitrogen bottom boundary.\n";

        return 1;
    }

    if (!approximately_equal(
            nitrogen(1, grid.ny - 1),
            1.0,
            1.0e-12
        ))
    {
        std::cerr
            << "Incorrect nitrogen top boundary.\n";

        return 1;
    }

    // ========================================================
    // Test 2: analytical advection
    // ========================================================

    counterflow::Field2D species_old(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D species_new(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D u(
        grid.nx,
        grid.ny,
        1.0
    );

    counterflow::Field2D v(
        grid.nx,
        grid.ny,
        0.0
    );

    // Y(x,y) = x
    //
    // dY/dx = 1
    // dY/dy = 0
    // Laplacian(Y) = 0
    //
    // With u = 1 and v = 0:
    //
    // Y_new = Y_old - dt

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            const std::size_t k =
                grid.index(i, j);

            species_old(i, j) =
                grid.x[k];
        }
    }

    const double dt = 0.1;
    const double diffusivity = 0.01;

    counterflow::advance_species_transport(
        species_old,
        u,
        v,
        species_new,
        grid,
        diffusivity,
        dt
    );

    // Center:
    //
    // Y_old = 0.5
    // Y_new = 0.5 - 0.1 = 0.4

    if (!approximately_equal(
            species_new(2, 2),
            0.4,
            1.0e-12
        ))
    {
        std::cerr
            << "Unexpected transported species value: "
            << species_new(2, 2)
            << '\n';

        return 1;
    }

    return 0;
}
