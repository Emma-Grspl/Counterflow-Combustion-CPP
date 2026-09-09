#include <cmath>
#include <iostream>
#include <stdexcept>

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
        5,
        3
    );

    // Grid dimensions
    if (grid.nx != 5 || grid.ny != 3)
    {
        std::cerr << "Unexpected grid dimensions.\n";
        return 1;
    }

    // Grid spacing
    if (!approximately_equal(grid.dx, 0.5, 1.0e-12))
    {
        std::cerr << "Unexpected dx value.\n";
        return 1;
    }

    if (!approximately_equal(grid.dy, 0.5, 1.0e-12))
    {
        std::cerr << "Unexpected dy value.\n";
        return 1;
    }

    // 2D -> 1D indexing
    if (grid.index(0, 0) != 0)
    {
        std::cerr << "Unexpected index for (0, 0).\n";
        return 1;
    }

    if (grid.index(2, 1) != 7)
    {
        std::cerr << "Unexpected index for (2, 1).\n";
        return 1;
    }

    if (grid.index(4, 2) != 14)
    {
        std::cerr << "Unexpected index for (4, 2).\n";
        return 1;
    }

    // Coordinates
    const std::size_t k = grid.index(2, 1);

    if (!approximately_equal(grid.x[k], 1.0, 1.0e-12))
    {
        std::cerr << "Unexpected x coordinate.\n";
        return 1;
    }

    if (!approximately_equal(grid.y[k], 0.5, 1.0e-12))
    {
        std::cerr << "Unexpected y coordinate.\n";
        return 1;
    }

    // Scalar field
    counterflow::Field2D field(5, 3, 0.0);

    if (field.nx() != 5 || field.ny() != 3)
    {
        std::cerr << "Unexpected field dimensions.\n";
        return 1;
    }

    field(2, 1) = 42.0;

    if (!approximately_equal(field(2, 1), 42.0, 1.0e-12))
    {
        std::cerr << "Field write/read failed.\n";
        return 1;
    }

    // Const access
    const counterflow::Field2D& const_field = field;

    if (!approximately_equal(const_field(2, 1), 42.0, 1.0e-12))
    {
        std::cerr << "Const field access failed.\n";
        return 1;
    }

    // Bounds checking
    bool caught_exception = false;

    try
    {
        [[maybe_unused]]
        const double value = field(5, 0);
    }
    catch (const std::out_of_range&)
    {
        caught_exception = true;
    }

    if (!caught_exception)
    {
        std::cerr << "Out-of-range access was not detected.\n";
        return 1;
    }

    return 0;
}
