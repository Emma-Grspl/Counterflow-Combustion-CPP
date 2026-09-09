#pragma once

#include <cstddef>

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include "counterflow/Grid.hpp"

namespace counterflow
{

void compute_pressure_rhs(
    const Field2D& u_star,
    const Field2D& v_star,
    Field2D& rhs,
    const Grid2D& grid,
    double rho,
    double dt
);


class PressurePoissonSolver
{
public:
    explicit PressurePoissonSolver(
        const Grid2D& grid
    );

    void solve(
        const Field2D& rhs,
        Field2D& pressure
    );

    [[nodiscard]]
    std::size_t system_size() const;

    [[nodiscard]]
    std::size_t nonzero_count() const;

private:
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using SparseSolver = Eigen::SparseLU<SparseMatrix>;

    std::size_t nx_;
    std::size_t ny_;

    double dx_;
    double dy_;

    SparseMatrix matrix_;
    SparseSolver solver_;

    Eigen::VectorXd rhs_vector_;
    Eigen::VectorXd solution_vector_;

    [[nodiscard]]
    std::size_t index(
        std::size_t i,
        std::size_t j
    ) const;

    void build_matrix();
};

} // namespace counterflow
