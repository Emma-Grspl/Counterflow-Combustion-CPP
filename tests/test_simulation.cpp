#include <cmath>
#include <iostream>

#include "counterflow/Config.hpp"
#include "counterflow/Simulation.hpp"

int main()
{
    counterflow::SimulationConfig config;

    // Small and fast integration test.
    config.nx = 7;
    config.ny = 7;

    config.final_time = 3.0e-8;
    config.nt = 4;

    counterflow::Simulation simulation(
        config
    );

    if (simulation.step_count() != 0)
    {
        std::cerr
            << "Simulation did not start at step zero.\n";

        return 1;
    }

    simulation.run(2);

    if (simulation.step_count() != 2)
    {
        std::cerr
            << "Incorrect simulation step count.\n";

        return 1;
    }

    const double expected_time =
        2.0 * config.dt();

    if (std::abs(
            simulation.time()
            - expected_time
        ) > 1.0e-15)
    {
        std::cerr
            << "Incorrect simulation time.\n";

        return 1;
    }

    if (!std::isfinite(
            simulation.max_temperature()
        ))
    {
        std::cerr
            << "Simulation produced a non-finite temperature.\n";

        return 1;
    }

    if (!(simulation.max_temperature() > 0.0))
    {
        std::cerr
            << "Simulation produced an invalid temperature.\n";

        return 1;
    }

    return 0;
}
