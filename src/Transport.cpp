#include "counterflow/Transport.hpp"

#include <algorithm>
#include <stdexcept>

namespace counterflow
{

void initialize_nitrogen_mass_fraction(
    Field2D& nitrogen,
    const Grid2D& grid
)
{
    if (nitrogen.nx() != grid.nx ||
        nitrogen.ny() != grid.ny)
    {
        throw std::invalid_argument(
            "Nitrogen field dimensions must match the grid."
        );
    }

    const std::size_t first_quarter =
        grid.nx / 4;

    const std::size_t half =
        grid.nx / 2;

    for (std::size_t i = 0; i < grid.nx; ++i)
    {
        if (i < first_quarter)
        {
            nitrogen(i, 0) = 0.79;
        }
        else if (i < half)
        {
            nitrogen(i, 0) = 1.0;
            nitrogen(i, grid.ny - 1) = 1.0;
        }
    }
}


void apply_nitrogen_boundary_conditions(
    Field2D& nitrogen,
    const Grid2D& grid
)
{
    if (nitrogen.nx() != grid.nx ||
        nitrogen.ny() != grid.ny)
    {
        throw std::invalid_argument(
            "Nitrogen field dimensions must match the grid."
        );
    }

    const std::size_t first_quarter =
        grid.nx / 4;

    const std::size_t half =
        grid.nx / 2;

    // Bottom and top: zero-gradient by default.
    for (std::size_t i = 0; i < grid.nx; ++i)
    {
        nitrogen(i, 0) =
            nitrogen(i, 1);

        nitrogen(i, grid.ny - 1) =
            nitrogen(i, grid.ny - 2);

        // Imposed inlet compositions.
        if (i < first_quarter)
        {
            nitrogen(i, 0) = 0.79;
        }
        else if (i < half)
        {
            nitrogen(i, 0) = 1.0;
            nitrogen(i, grid.ny - 1) = 1.0;
        }
    }

    // Left and right: zero-gradient.
    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        nitrogen(0, j) =
            nitrogen(1, j);

        nitrogen(grid.nx - 1, j) =
            nitrogen(grid.nx - 2, j);
    }
}


void advance_species_transport(
    const Field2D& species_old,
    const Field2D& u,
    const Field2D& v,
    Field2D& species_new,
    const Grid2D& grid,
    double diffusivity,
    double dt
)
{
    const bool dimensions_match =
        species_old.nx() == grid.nx &&
        species_old.ny() == grid.ny &&
        species_new.nx() == grid.nx &&
        species_new.ny() == grid.ny &&
        u.nx() == grid.nx &&
        u.ny() == grid.ny &&
        v.nx() == grid.nx &&
        v.ny() == grid.ny;

    if (!dimensions_match)
    {
        throw std::invalid_argument(
            "Species transport fields must match the grid."
        );
    }

    if (diffusivity < 0.0)
    {
        throw std::invalid_argument(
            "Diffusivity cannot be negative."
        );
    }

    if (dt <= 0.0)
    {
        throw std::invalid_argument(
            "Time step must be strictly positive."
        );
    }

    const double inv_2dx =
        1.0 / (2.0 * grid.dx);

    const double inv_2dy =
        1.0 / (2.0 * grid.dy);

    const double inv_dx2 =
        1.0 / (grid.dx * grid.dx);

    const double inv_dy2 =
        1.0 / (grid.dy * grid.dy);

    // Preserve the current field first.
    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            species_new(i, j) =
                species_old(i, j);
        }
    }

    // Explicit centered advection-diffusion update.
    for (std::size_t j = 1; j < grid.ny - 1; ++j)
    {
        for (std::size_t i = 1; i < grid.nx - 1; ++i)
        {
            const double dY_dx =
                (species_old(i + 1, j)
                 - species_old(i - 1, j))
                * inv_2dx;

            const double dY_dy =
                (species_old(i, j + 1)
                 - species_old(i, j - 1))
                * inv_2dy;

            const double laplacian =
                (species_old(i + 1, j)
                 - 2.0 * species_old(i, j)
                 + species_old(i - 1, j))
                * inv_dx2
                +
                (species_old(i, j + 1)
                 - 2.0 * species_old(i, j)
                 + species_old(i, j - 1))
                * inv_dy2;

            const double advection =
                u(i, j) * dY_dx
                + v(i, j) * dY_dy;

            species_new(i, j) =
                species_old(i, j)
                - dt * advection
                + diffusivity * dt * laplacian;
        }
    }

    // Physical mass fractions must remain in [0, 1].
    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            species_new(i, j) =
                std::clamp(
                    species_new(i, j),
                    0.0,
                    1.0
                );
        }
    }
}

} // namespace counterflow
