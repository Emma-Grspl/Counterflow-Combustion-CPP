#include "counterflow/BoundaryConditions.hpp"

#include <stdexcept>

namespace counterflow
{

void initialize_vertical_velocity(
    Field2D& v,
    const Grid2D& grid
)
{
    if (v.nx() != grid.nx || v.ny() != grid.ny)
    {
        throw std::invalid_argument(
            "Velocity field dimensions must match the grid."
        );
    }

    const std::size_t fast_jet_end = grid.nx / 4;
    const std::size_t slow_jet_end = grid.nx / 2;

    for (std::size_t i = 0; i < grid.nx; ++i)
    {
        // Bottom boundary
        if (i < fast_jet_end)
        {
            v(i, 0) = 1.0;
        }
        else if (i < slow_jet_end)
        {
            v(i, 0) = 0.2;
        }

        // Top boundary
        if (i < fast_jet_end)
        {
            v(i, grid.ny - 1) = -1.0;
        }
        else if (i < slow_jet_end)
        {
            v(i, grid.ny - 1) = -0.2;
        }
    }
}

} // namespace counterflow
