#pragma once

#include "counterflow/Grid.hpp"

namespace counterflow
{

void initialize_vertical_velocity(
    Field2D& v,
    const Grid2D& grid
);

} // namespace counterflow
