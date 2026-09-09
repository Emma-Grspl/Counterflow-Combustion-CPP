#include "counterflow/PressurePoisson.hpp"

#include <stdexcept>

namespace counterflow
{

void compute_pressure_rhs(
    const Field2D& u_star,
    const Field2D& v_star,
    Field2D& rhs,
    const Grid2D& grid,
    double rho,
    double dt
)
{
    const bool dimensions_match =
        u_star.nx() == grid.nx &&
        u_star.ny() == grid.ny &&
        v_star.nx() == grid.nx &&
        v_star.ny() == grid.ny &&
        rhs.nx() == grid.nx &&
        rhs.ny() == grid.ny;

    if (!dimensions_match)
    {
        throw std::invalid_argument(
            "Pressure RHS field dimensions must match the grid."
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

    const double pressure_scale =
        rho / dt;

    for (std::size_t j = 1; j < grid.ny - 1; ++j)
    {
        for (std::size_t i = 1; i < grid.nx - 1; ++i)
        {
            const double du_dx =
                (u_star(i + 1, j)
                 - u_star(i - 1, j))
                * inv_2dx;

            const double dv_dy =
                (v_star(i, j + 1)
                 - v_star(i, j - 1))
                * inv_2dy;

            rhs(i, j) =
                pressure_scale
                * (du_dx + dv_dy);
        }
    }
}

} // namespace counterflow
