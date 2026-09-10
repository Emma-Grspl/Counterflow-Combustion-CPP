#pragma once

#include <cstddef>

#include "counterflow/Combustion.hpp"
#include "counterflow/Grid.hpp"

namespace counterflow
{

void initialize_reactive_fields(
    Field2D& ch4,
    Field2D& o2,
    Field2D& h2o,
    Field2D& co2,
    Field2D& temperature,
    const Grid2D& grid
);

void update_nitrogen_from_mass_closure(
    Field2D& nitrogen,
    const Field2D& ch4,
    const Field2D& o2,
    const Field2D& h2o,
    const Field2D& co2,
    const Grid2D& grid
);

void apply_reactive_boundary_conditions(
    Field2D& ch4,
    Field2D& o2,
    Field2D& h2o,
    Field2D& co2,
    Field2D& temperature,
    const Grid2D& grid
);


class ReactiveTransportStepper
{
public:
    ReactiveTransportStepper(
        const Grid2D& grid,
        double density,
        double heat_capacity,
        double diffusivity,
        double hydro_dt,
        double maximum_temperature
    );

    void advance(
        Field2D& ch4,
        Field2D& o2,
        Field2D& h2o,
        Field2D& co2,
        Field2D& temperature,
        const Field2D& u,
        const Field2D& v,
        const Grid2D& grid,
        bool evolve_temperature = true
    );

    [[nodiscard]]
    std::size_t substeps() const;

    [[nodiscard]]
    double chemical_dt() const;

private:
    std::size_t nx_;
    std::size_t ny_;

    double density_;
    double heat_capacity_;
    double diffusivity_;
    double hydro_dt_;

    ChemicalSubcycling subcycling_;

    Field2D ch4_buffer_;
    Field2D o2_buffer_;
    Field2D h2o_buffer_;
    Field2D co2_buffer_;
    Field2D temperature_buffer_;
};

} // namespace counterflow
