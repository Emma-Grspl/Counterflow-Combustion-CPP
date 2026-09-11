#include "counterflow/BoundaryConditions.hpp"

#include <stdexcept>


namespace counterflow
{

namespace
{

double imposed_vertical_velocity(
    std::size_t i,
    std::size_t nx,
    bool top
)
{
    const std::size_t fast_jet_end =
        nx / 4;

    const std::size_t slow_jet_end =
        nx / 2;

    double magnitude = 0.0;

    if (i < fast_jet_end)
    {
        magnitude = 1.0;
    }
    else if (i < slow_jet_end)
    {
        magnitude = 0.2;
    }

    return top
        ? -magnitude
        : magnitude;
}


void validate_velocity_field(
    const Field2D& field,
    const Grid2D& grid
)
{
    if (
        field.nx() != grid.nx
        || field.ny() != grid.ny
    )
    {
        throw std::invalid_argument(
            "Velocity field dimensions must match the grid."
        );
    }
}

} // namespace


void initialize_vertical_velocity(
    Field2D& v,
    const Grid2D& grid
)
{
    validate_velocity_field(
        v,
        grid
    );

    for (
        std::size_t i = 0;
        i < grid.nx;
        ++i
    )
    {
        v(i, 0) =
            imposed_vertical_velocity(
                i,
                grid.nx,
                false
            );

        v(i, grid.ny - 1) =
            imposed_vertical_velocity(
                i,
                grid.nx,
                true
            );
    }
}


void apply_velocity_boundary_conditions(
    Field2D& u,
    Field2D& v,
    const Grid2D& grid
)
{
    validate_velocity_field(
        u,
        grid
    );

    validate_velocity_field(
        v,
        grid
    );

    // --------------------------------------------------------
    // Left boundary: slip wall
    //
    // u = 0
    // dv/dx = 0
    // --------------------------------------------------------

    for (
        std::size_t j = 1;
        j < grid.ny - 1;
        ++j
    )
    {
        u(0, j) = 0.0;

        v(0, j) =
            v(1, j);
    }

    // --------------------------------------------------------
    // Right boundary: open outlet
    //
    // du/dx = 0
    // dv/dx = 0
    //
    // Pressure is imposed separately as p = 0 by the
    // pressure Poisson solver.
    // --------------------------------------------------------

    const std::size_t right =
        grid.nx - 1;

    for (
        std::size_t j = 1;
        j < grid.ny - 1;
        ++j
    )
    {
        u(right, j) =
            u(right - 1, j);

        v(right, j) =
            v(right - 1, j);
    }

    // --------------------------------------------------------
    // Bottom / top boundaries
    //
    // Horizontal velocity is zero.
    // Vertical velocity contains the imposed counterflow jets.
    //
    // These are applied last so the injection conditions own
    // the corner values.
    // --------------------------------------------------------

    const std::size_t top =
        grid.ny - 1;

    for (
        std::size_t i = 0;
        i < grid.nx;
        ++i
    )
    {
        u(i, 0) = 0.0;
        u(i, top) = 0.0;

        v(i, 0) =
            imposed_vertical_velocity(
                i,
                grid.nx,
                false
            );

        v(i, top) =
            imposed_vertical_velocity(
                i,
                grid.nx,
                true
            );
    }
}

} // namespace counterflow
