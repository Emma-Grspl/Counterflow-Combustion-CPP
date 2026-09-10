#pragma once

#include <cstddef>
#include <optional>

namespace counterflow
{

struct SimulationConfig
{
    // Physical parameters
    double rho = 1.1614;      // kg/m^3
    double nu = 15.0e-6;      // m^2/s
    double diffusivity = 15.0e-6; // m^2/s
    double cp = 1200.0;       // J/(kg.K)

    // Domain dimensions
    double lx = 0.002;        // m
    double ly = 0.002;        // m

    // Grid resolution
    std::size_t nx = 50;
    std::size_t ny = 50;

    // Time integration
    double final_time = 0.01; // s
    std::size_t nt = 3000;

    // Flow steady-state detection
    double steady_state_tolerance = 1.0e-8;

    // Reference reproduction:
    // if a value is present, activate the energy equation
    // after this prescribed hydrodynamic step.
    //
    // std::nullopt selects automatic steady-state detection.
    std::optional<std::size_t>
        prescribed_energy_activation_step = 1524;

    [[nodiscard]]
    double dt() const
    {
        return final_time / static_cast<double>(nt - 1);
    }
};

} // namespace counterflow
