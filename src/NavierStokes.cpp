#include "counterflow/NavierStokes.hpp"

#include <stdexcept>

namespace counterflow
{

void compute_intermediate_velocity(
    const Field2D& u_old,
    const Field2D& v_old,
    Field2D& u_star,
    Field2D& v_star,
    const Grid2D& grid,
    double dt,
    double nu
)
{
    const bool dimensions_match =
        u_old.nx() == grid.nx &&
        u_old.ny() == grid.ny &&
        v_old.nx() == grid.nx &&
        v_old.ny() == grid.ny &&
        u_star.nx() == grid.nx &&
        u_star.ny() == grid.ny &&
        v_star.nx() == grid.nx &&
        v_star.ny() == grid.ny;

    if (!dimensions_match)
    {
        throw std::invalid_argument(
            "Velocity field dimensions must match the grid."
        );
    }

    const double inv_2dx = 1.0 / (2.0 * grid.dx);
    const double inv_2dy = 1.0 / (2.0 * grid.dy);

    const double inv_dx2 = 1.0 / (grid.dx * grid.dx);
    const double inv_dy2 = 1.0 / (grid.dy * grid.dy);

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            const bool boundary =
                i == 0 ||
                i == grid.nx - 1 ||
                j == 0 ||
                j == grid.ny - 1;

            if (boundary)
            {
                u_star(i, j) = u_old(i, j);
                v_star(i, j) = v_old(i, j);
                continue;
            }

            // First-order spatial derivatives
            const double du_dx =
                (u_old(i + 1, j) - u_old(i - 1, j))
                * inv_2dx;

            const double du_dy =
                (u_old(i, j + 1) - u_old(i, j - 1))
                * inv_2dy;

            const double dv_dx =
                (v_old(i + 1, j) - v_old(i - 1, j))
                * inv_2dx;

            const double dv_dy =
                (v_old(i, j + 1) - v_old(i, j - 1))
                * inv_2dy;

            // Laplacians
            const double laplacian_u =
                (u_old(i + 1, j)
                 - 2.0 * u_old(i, j)
                 + u_old(i - 1, j))
                * inv_dx2
                +
                (u_old(i, j + 1)
                 - 2.0 * u_old(i, j)
                 + u_old(i, j - 1))
                * inv_dy2;

            const double laplacian_v =
                (v_old(i + 1, j)
                 - 2.0 * v_old(i, j)
                 + v_old(i - 1, j))
                * inv_dx2
                +
                (v_old(i, j + 1)
                 - 2.0 * v_old(i, j)
                 + v_old(i, j - 1))
                * inv_dy2;

            // Convective terms
            const double advection_u =
                u_old(i, j) * du_dx
                + v_old(i, j) * du_dy;

            const double advection_v =
                u_old(i, j) * dv_dx
                + v_old(i, j) * dv_dy;

            // Predictor step
            u_star(i, j) =
                u_old(i, j)
                - dt * advection_u
                + nu * dt * laplacian_u;

            v_star(i, j) =
                v_old(i, j)
                - dt * advection_v
                + nu * dt * laplacian_v;
        }
    }
}

} // namespace counterflow
