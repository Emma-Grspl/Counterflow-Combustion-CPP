#include <cmath>
#include <iostream>

#include "counterflow/Config.hpp"

bool approximately_equal(double a, double b, double tolerance)
{
    return std::abs(a - b) < tolerance;
}

int main()
{
    const counterflow::SimulationConfig config;

    if (config.nx != 50)
    {
        std::cerr << "Expected nx = 50, got " << config.nx << '\n';
        return 1;
    }

    if (config.ny != 50)
    {
        std::cerr << "Expected ny = 50, got " << config.ny << '\n';
        return 1;
    }

    if (config.nt != 3000)
    {
        std::cerr << "Expected nt = 3000, got " << config.nt << '\n';
        return 1;
    }

    if (!approximately_equal(config.rho, 1.1614, 1.0e-12))
    {
        std::cerr << "Unexpected rho value\n";
        return 1;
    }

    if (!approximately_equal(config.nu, 15.0e-6, 1.0e-12))
    {
        std::cerr << "Unexpected nu value\n";
        return 1;
    }

    if (!approximately_equal(config.diffusivity, 15.0e-6, 1.0e-12))
    {
        std::cerr << "Unexpected diffusivity value\n";
        return 1;
    }

    const double expected_dt =
        config.final_time / static_cast<double>(config.nt - 1);

    if (!approximately_equal(config.dt(), expected_dt, 1.0e-15))
    {
        std::cerr << "Unexpected time-step value\n";
        return 1;
    }

    return 0;
}
