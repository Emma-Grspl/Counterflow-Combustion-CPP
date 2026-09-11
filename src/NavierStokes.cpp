#include "counterflow/NavierStokes.hpp"
#include "counterflow/BoundaryConditions.hpp"

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


void correct_velocity(
    const Field2D& u_star,
    const Field2D& v_star,
    const Field2D& pressure,
    Field2D& u_new,
    Field2D& v_new,
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
        pressure.nx() == grid.nx &&
        pressure.ny() == grid.ny &&
        u_new.nx() == grid.nx &&
        u_new.ny() == grid.ny &&
        v_new.nx() == grid.nx &&
        v_new.ny() == grid.ny;

    if (!dimensions_match)
    {
        throw std::invalid_argument(
            "Velocity correction fields must match the grid."
        );
    }

    if (rho <= 0.0)
    {
        throw std::invalid_argument(
            "Density must be strictly positive."
        );
    }

    if (dt <= 0.0)
    {
        throw std::invalid_argument(
            "Time step must be strictly positive."
        );
    }

    const double pressure_scale =
        dt / rho;

    const double inv_dx =
        1.0 / grid.dx;

    const double inv_dy =
        1.0 / grid.dy;

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            if (i == 0 ||
                i == grid.nx - 1 ||
                j == 0 ||
                j == grid.ny - 1)
            {
                u_new(i, j) =
                    u_star(i, j);

                v_new(i, j) =
                    v_star(i, j);

                continue;
            }

            // Forward pressure gradient G^+.
            const double dp_dx =
                (pressure(i + 1, j)
                 - pressure(i, j))
                * inv_dx;

            const double dp_dy =
                (pressure(i, j + 1)
                 - pressure(i, j))
                * inv_dy;

            u_new(i, j) =
                u_star(i, j)
                - pressure_scale * dp_dx;

            v_new(i, j) =
                v_star(i, j)
                - pressure_scale * dp_dy;
        }
    }
}


NavierStokesStepper::NavierStokesStepper(
    const Grid2D& grid,
    double rho,
    double nu,
    double dt
)
    : nx_(grid.nx),
      ny_(grid.ny),
      rho_(rho),
      nu_(nu),
      dt_(dt),
      u_star_(grid.nx, grid.ny, 0.0),
      v_star_(grid.nx, grid.ny, 0.0),
      pressure_rhs_(grid.nx, grid.ny, 0.0),
      pressure_solver_(grid)
{
    if (rho_ <= 0.0)
    {
        throw std::invalid_argument(
            "Density must be strictly positive."
        );
    }

    if (nu_ < 0.0)
    {
        throw std::invalid_argument(
            "Kinematic viscosity cannot be negative."
        );
    }

    if (dt_ <= 0.0)
    {
        throw std::invalid_argument(
            "Time step must be strictly positive."
        );
    }
}


void NavierStokesStepper::advance(
    const Field2D& u_old,
    const Field2D& v_old,
    Field2D& u_new,
    Field2D& v_new,
    Field2D& pressure,
    const Grid2D& grid
)
{
    if (grid.nx != nx_ || grid.ny != ny_)
    {
        throw std::invalid_argument(
            "Grid dimensions do not match the Navier-Stokes stepper."
        );
    }

    compute_intermediate_velocity(
        u_old,
        v_old,
        u_star_,
        v_star_,
        grid,
        dt_,
        nu_
    );

    apply_velocity_boundary_conditions(
        u_star_,
        v_star_,
        grid
    );

    compute_pressure_rhs(
        u_star_,
        v_star_,
        pressure_rhs_,
        grid,
        rho_,
        dt_
    );

    pressure_solver_.solve(
        pressure_rhs_,
        pressure
    );

    correct_velocity(
        u_star_,
        v_star_,
        pressure,
        u_new,
        v_new,
        grid,
        rho_,
        dt_
    );

    apply_velocity_boundary_conditions(
        u_new,
        v_new,
        grid
    );
}

} // namespace counterflow
