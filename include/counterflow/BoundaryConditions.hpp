#pragma once

#include "counterflow/Grid.hpp"

namespace counterflow
{

void initialize_vertical_velocity(
    Field2D& v,
    const Grid2D& grid
);

void apply_velocity_boundary_conditions(
    Field2D& u,
    Field2D& v,
    const Grid2D& grid
);

} // namespace counterflow
