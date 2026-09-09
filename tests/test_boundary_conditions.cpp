#include <cmath>
#include <iostream>

#include "counterflow/BoundaryConditions.hpp"
#include "counterflow/Grid.hpp"

bool approximately_equal(double a, double b, double tolerance)
{
    return std::abs(a - b) < tolerance;
}

int main()
{
    const counterflow::Grid2D grid(
        2.0,
        1.0,
        8,
        5
    );

    counterflow::Field2D v(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::initialize_vertical_velocity(v, grid);

    // Bottom boundary
    if (!approximately_equal(v(0, 0), 1.0, 1.0e-12) ||
        !approximately_equal(v(1, 0), 1.0, 1.0e-12) ||
        !approximately_equal(v(2, 0), 0.2, 1.0e-12) ||
        !approximately_equal(v(3, 0), 0.2, 1.0e-12))
    {
        std::cerr << "Incorrect bottom velocity boundary.\n";
        return 1;
    }

    // Top boundary
    const std::size_t top = grid.ny - 1;

    if (!approximately_equal(v(0, top), -1.0, 1.0e-12) ||
        !approximately_equal(v(1, top), -1.0, 1.0e-12) ||
        !approximately_equal(v(2, top), -0.2, 1.0e-12) ||
        !approximately_equal(v(3, top), -0.2, 1.0e-12))
    {
        std::cerr << "Incorrect top velocity boundary.\n";
        return 1;
    }

    // Interior must remain zero
    if (!approximately_equal(v(3, 2), 0.0, 1.0e-12))
    {
        std::cerr << "Interior velocity was unexpectedly modified.\n";
        return 1;
    }

    return 0;
}
