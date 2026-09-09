#pragma once

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

} // namespace counterflow
