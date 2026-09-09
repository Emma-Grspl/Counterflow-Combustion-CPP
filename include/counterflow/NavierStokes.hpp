#pragma once

#include "counterflow/Grid.hpp"

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

} // namespace counterflow
