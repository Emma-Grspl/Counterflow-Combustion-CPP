#include <algorithm>
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


    // ========================================================
    // Mixture mass-fraction closure
    // ========================================================

    const counterflow::Grid2D& grid =
        simulation.grid();

    const counterflow::Field2D& nitrogen =
        simulation.nitrogen();

    const counterflow::Field2D& ch4 =
        simulation.ch4();

    const counterflow::Field2D& o2 =
        simulation.o2();

    const counterflow::Field2D& h2o =
        simulation.h2o();

    const counterflow::Field2D& co2 =
        simulation.co2();

    double maximum_mass_error = 0.0;

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            const double mass_sum =
                nitrogen(i, j)
                + ch4(i, j)
                + o2(i, j)
                + h2o(i, j)
                + co2(i, j);

            maximum_mass_error =
                std::max(
                    maximum_mass_error,
                    std::abs(mass_sum - 1.0)
                );
        }
    }

    if (maximum_mass_error > 1.0e-12)
    {
        std::cerr
            << "Mass-fraction closure violated: "
            << maximum_mass_error
            << '\n';

        return 1;
    }

    return 0;
}
