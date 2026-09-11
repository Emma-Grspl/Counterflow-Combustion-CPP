#include <charconv>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <string_view>

#include "counterflow/Config.hpp"
#include "counterflow/Output.hpp"
#include "counterflow/Simulation.hpp"


namespace
{

bool parse_size(
    std::string_view text,
    std::size_t& value
)
{
    const char* begin = text.data();
    const char* end = begin + text.size();

    const auto result = std::from_chars(
        begin,
        end,
        value
    );

    return (
        result.ec == std::errc{}
        && result.ptr == end
    );
}

} // namespace


int main(
    int argc,
    char* argv[]
)
{
    if (argc != 4)
    {
        std::cerr
            << "Usage: "
            << argv[0]
            << " <grid_points> <time_levels> <output.csv>\n";

        return 1;
    }

    std::size_t grid_points = 0;
    std::size_t time_levels = 0;

    if (
        !parse_size(argv[1], grid_points)
        || grid_points < 3
    )
    {
        std::cerr
            << "Invalid grid size: "
            << argv[1]
            << '\n';

        return 1;
    }

    if (
        !parse_size(argv[2], time_levels)
        || time_levels < 2
    )
    {
        std::cerr
            << "Invalid number of time levels: "
            << argv[2]
            << '\n';

        return 1;
    }

    counterflow::SimulationConfig config;

    config.nx = grid_points;
    config.ny = grid_points;
    config.nt = time_levels;

    // Physical ignition/energy-activation time from the
    // original reference calculation:
    //
    //     step = 1524
    //     dt   = 0.01 / 2999
    //
    // When dt changes, preserve this physical time rather
    // than preserving the integer step number.
    constexpr double reference_final_time = 0.01;
    constexpr std::size_t reference_intervals = 2999;
    constexpr std::size_t reference_activation_step = 1524;

    const double reference_dt =
        reference_final_time
        / static_cast<double>(
            reference_intervals
        );

    const double activation_time =
        static_cast<double>(
            reference_activation_step
        )
        * reference_dt;

    const auto activation_step =
        static_cast<std::size_t>(
            std::llround(
                activation_time
                / config.dt()
            )
        );

    config.prescribed_energy_activation_step =
        activation_step;

    const std::filesystem::path output_path{
        argv[3]
    };

    std::cout
        << "Counterflow convergence run\n"
        << "===========================\n"
        << "Grid              : "
        << config.nx
        << " x "
        << config.ny
        << '\n'
        << "Time levels       : "
        << config.nt
        << '\n'
        << "Time steps        : "
        << config.nt - 1
        << '\n'
        << "dt                : "
        << config.dt()
        << " s\n"
        << "Activation step   : "
        << activation_step
        << '\n'
        << "Activation time   : "
        << activation_step * config.dt()
        << " s\n\n";

    counterflow::Simulation simulation{
        config
    };

    simulation.run(
        config.nt - 1
    );

    counterflow::write_simulation_csv(
        simulation,
        output_path
    );

    std::cout
        << "Final Tmax        : "
        << simulation.max_temperature()
        << " K\n"
        << "Output            : "
        << output_path
        << '\n';

    return 0;
}
