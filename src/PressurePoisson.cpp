#include "counterflow/PressurePoisson.hpp"

#include <stdexcept>
#include <vector>

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
            "Pressure RHS fields must match the grid."
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

    // Boundary rows of the Poisson system represent
    // homogeneous pressure boundary conditions.
    // Reset the full RHS so stale values cannot survive
    // from a previous time step.
    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            rhs(i, j) = 0.0;
        }
    }

    const double inv_dx =
        1.0 / grid.dx;

    const double inv_dy =
        1.0 / grid.dy;

    const double pressure_scale =
        rho / dt;

    // Backward divergence D^-.
    //
    // Combined with the forward pressure gradient G^+
    // used in correct_velocity(), D^- G^+ gives exactly
    // the standard five-point Laplacian.
    for (std::size_t j = 1; j < grid.ny - 1; ++j)
    {
        for (std::size_t i = 1; i < grid.nx - 1; ++i)
        {
            const double du_dx =
                (u_star(i, j)
                 - u_star(i - 1, j))
                * inv_dx;

            const double dv_dy =
                (v_star(i, j)
                 - v_star(i, j - 1))
                * inv_dy;

            rhs(i, j) =
                pressure_scale
                * (du_dx + dv_dy);
        }
    }
}


PressurePoissonSolver::PressurePoissonSolver(
    const Grid2D& grid
)
    : nx_(grid.nx),
      ny_(grid.ny),
      dx_(grid.dx),
      dy_(grid.dy),
      matrix_(
          static_cast<Eigen::Index>(grid.nx * grid.ny),
          static_cast<Eigen::Index>(grid.nx * grid.ny)
      ),
      rhs_vector_(
          static_cast<Eigen::Index>(grid.nx * grid.ny)
      ),
      solution_vector_(
          static_cast<Eigen::Index>(grid.nx * grid.ny)
      )
{
    build_matrix();

    solver_.analyzePattern(matrix_);

    solver_.factorize(matrix_);

    if (solver_.info() != Eigen::Success)
    {
        throw std::runtime_error(
            "Pressure Poisson matrix factorization failed."
        );
    }
}


std::size_t PressurePoissonSolver::index(
    std::size_t i,
    std::size_t j
) const
{
    return i + j * nx_;
}


void PressurePoissonSolver::build_matrix()
{
    using Triplet = Eigen::Triplet<double>;

    const std::size_t n =
        nx_ * ny_;

    std::vector<Triplet> triplets;

    triplets.reserve(5 * n);

    const double inv_dx =
        1.0 / dx_;

    const double inv_dy =
        1.0 / dy_;

    const double inv_dx2 =
        1.0 / (dx_ * dx_);

    const double inv_dy2 =
        1.0 / (dy_ * dy_);

    for (std::size_t j = 0; j < ny_; ++j)
    {
        for (std::size_t i = 0; i < nx_; ++i)
        {
            const Eigen::Index row =
                static_cast<Eigen::Index>(
                    index(i, j)
                );

            // Left boundary:
            // dp/dx = 0
            if (i == 0)
            {
                triplets.emplace_back(
                    row,
                    row,
                    -inv_dx
                );

                triplets.emplace_back(
                    row,
                    static_cast<Eigen::Index>(
                        index(i + 1, j)
                    ),
                    inv_dx
                );

                continue;
            }

            // Right boundary:
            // p = 0
            if (i == nx_ - 1)
            {
                triplets.emplace_back(
                    row,
                    row,
                    1.0
                );

                continue;
            }

            // Bottom boundary:
            // dp/dy = 0
            if (j == 0)
            {
                triplets.emplace_back(
                    row,
                    row,
                    -inv_dy
                );

                triplets.emplace_back(
                    row,
                    static_cast<Eigen::Index>(
                        index(i, j + 1)
                    ),
                    inv_dy
                );

                continue;
            }

            // Top boundary:
            // dp/dy = 0
            if (j == ny_ - 1)
            {
                triplets.emplace_back(
                    row,
                    row,
                    inv_dy
                );

                triplets.emplace_back(
                    row,
                    static_cast<Eigen::Index>(
                        index(i, j - 1)
                    ),
                    -inv_dy
                );

                continue;
            }

            // Interior 5-point Laplacian
            triplets.emplace_back(
                row,
                row,
                -2.0 * inv_dx2
                - 2.0 * inv_dy2
            );

            triplets.emplace_back(
                row,
                static_cast<Eigen::Index>(
                    index(i + 1, j)
                ),
                inv_dx2
            );

            triplets.emplace_back(
                row,
                static_cast<Eigen::Index>(
                    index(i - 1, j)
                ),
                inv_dx2
            );

            triplets.emplace_back(
                row,
                static_cast<Eigen::Index>(
                    index(i, j + 1)
                ),
                inv_dy2
            );

            triplets.emplace_back(
                row,
                static_cast<Eigen::Index>(
                    index(i, j - 1)
                ),
                inv_dy2
            );
        }
    }

    matrix_.setFromTriplets(
        triplets.begin(),
        triplets.end()
    );

    matrix_.makeCompressed();
}


void PressurePoissonSolver::solve(
    const Field2D& rhs,
    Field2D& pressure
)
{
    const bool dimensions_match =
        rhs.nx() == nx_ &&
        rhs.ny() == ny_ &&
        pressure.nx() == nx_ &&
        pressure.ny() == ny_;

    if (!dimensions_match)
    {
        throw std::invalid_argument(
            "Pressure solve field dimensions do not match."
        );
    }

    for (std::size_t j = 0; j < ny_; ++j)
    {
        for (std::size_t i = 0; i < nx_; ++i)
        {
            const Eigen::Index k =
                static_cast<Eigen::Index>(
                    index(i, j)
                );

            rhs_vector_[k] =
                rhs(i, j);
        }
    }

    solution_vector_ =
        solver_.solve(rhs_vector_);

    if (solver_.info() != Eigen::Success)
    {
        throw std::runtime_error(
            "Pressure Poisson solve failed."
        );
    }

    for (std::size_t j = 0; j < ny_; ++j)
    {
        for (std::size_t i = 0; i < nx_; ++i)
        {
            const Eigen::Index k =
                static_cast<Eigen::Index>(
                    index(i, j)
                );

            pressure(i, j) =
                solution_vector_[k];
        }
    }
}


std::size_t PressurePoissonSolver::system_size() const
{
    return nx_ * ny_;
}


std::size_t PressurePoissonSolver::nonzero_count() const
{
    return static_cast<std::size_t>(
        matrix_.nonZeros()
    );
}

} // namespace counterflow
