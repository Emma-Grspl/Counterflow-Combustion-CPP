#pragma once

#include "counterflow/Grid.hpp"

namespace counterflow
{

void initialize_nitrogen_mass_fraction(
    Field2D& nitrogen,
    const Grid2D& grid
);

void apply_nitrogen_boundary_conditions(
    Field2D& nitrogen,
    const Grid2D& grid
);

void advance_scalar_transport(
    const Field2D& scalar_old,
    const Field2D& u,
    const Field2D& v,
    Field2D& scalar_new,
    const Grid2D& grid,
    double diffusivity,
    double dt
);

void advance_species_transport(
    const Field2D& species_old,
    const Field2D& u,
    const Field2D& v,
    Field2D& species_new,
    const Grid2D& grid,
    double diffusivity,
    double dt
);

} // namespace counterflow
