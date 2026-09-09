#pragma once

#include <cstddef>

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

    [[nodiscard]]
    double dt() const
    {
        return final_time / static_cast<double>(nt - 1);
    }
};

} // namespace counterflow
