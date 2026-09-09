#pragma once

#include "counterflow/Grid.hpp"
#include "counterflow/PressurePoisson.hpp"

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
);

void correct_velocity(
    const Field2D& u_star,
    const Field2D& v_star,
    const Field2D& pressure,
    Field2D& u_new,
    Field2D& v_new,
    const Grid2D& grid,
    double rho,
    double dt
);


class NavierStokesStepper
{
public:
    NavierStokesStepper(
        const Grid2D& grid,
        double rho,
        double nu,
        double dt
    );

    void advance(
        const Field2D& u_old,
        const Field2D& v_old,
        Field2D& u_new,
        Field2D& v_new,
        Field2D& pressure,
        const Grid2D& grid
    );

private:
    std::size_t nx_;
    std::size_t ny_;

    double rho_;
    double nu_;
    double dt_;

    Field2D u_star_;
    Field2D v_star_;
    Field2D pressure_rhs_;

    PressurePoissonSolver pressure_solver_;
};

} // namespace counterflow
