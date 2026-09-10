#include <charconv>
#include <cstddef>
#include <iostream>
#include <optional>
#include <filesystem>
#include <string_view>

#include "counterflow/Config.hpp"
#include "counterflow/Simulation.hpp"
#include "counterflow/Output.hpp"

namespace
{

bool parse_step_count(
    std::string_view text,
    std::size_t& value
)
{
    const char* begin =
        text.data();

    const char* end =
        text.data() + text.size();

    const auto result =
        std::from_chars(
            begin,
            end,
            value
        );

    return result.ec == std::errc{}
        && result.ptr == end;
}

} // namespace


int main(
    int argc,
    char* argv[]
)
{
    counterflow::SimulationConfig config;

    std::size_t number_of_steps = 10;

    if (argc > 4)
    {
        std::cerr
            << "Usage: "
            << argv[0]
            << " [number_of_steps|all] "
               "[output.csv] [reference|auto]\n";

        return 1;
    }

    if (argc >= 2)
    {
        const std::string_view argument{
            argv[1]
        };

        if (argument == "all")
        {
            number_of_steps =
                config.nt - 1;
        }
        else if (!parse_step_count(
                     argument,
                     number_of_steps
                 ))
        {
            std::cerr
                << "Invalid number of steps: "
                << argument
                << '\n';

            return 1;
        }
    }

    if (argc == 4)
    {
        const std::string_view mode{
            argv[3]
        };

        if (mode == "reference")
        {
            config.prescribed_energy_activation_step =
                1524;
        }
        else if (mode == "auto")
        {
            config.prescribed_energy_activation_step =
                std::nullopt;
        }
        else
        {
            std::cerr
                << "Invalid activation mode: "
                << mode
                << "\nExpected 'reference' or 'auto'.\n";

            return 1;
        }
    }

    if (number_of_steps > config.nt - 1)
    {
        std::cerr
            << "Requested "
            << number_of_steps
            << " steps, but the configuration "
               "allows at most "
            << config.nt - 1
            << ".\n";

        return 1;
    }

    std::cout
        << "Counterflow Combustion C++ Solver\n"
        << "---------------------------------\n"
        << "Grid: "
        << config.nx
        << " x "
        << config.ny
        << '\n'
        << "dt: "
        << config.dt()
        << " s\n"
        << "Requested steps: "
        << number_of_steps
        << '\n'
        << "Energy activation: "
        << (
            config.prescribed_energy_activation_step.has_value()
                ? "reference step"
                : "automatic steady state"
           )
        << "\n\n";

    counterflow::Simulation simulation(
        config
    );

    simulation.run(
        number_of_steps
    );

    if (argc == 3)
    {
        const std::filesystem::path output_path{
            argv[2]
        };

        counterflow::write_simulation_csv(
            simulation,
            output_path
        );

        std::cout
            << "Output written to: "
            << output_path
            << '\n';
    }

    std::cout
        << "Simulation completed.\n"
        << "Steps: "
        << simulation.step_count()
        << '\n'
        << "Physical time: "
        << simulation.time()
        << " s\n"
        << "Maximum temperature: "
        << simulation.max_temperature()
        << " K\n";

    if (simulation.steady_state_detected())
    {
        std::cout
            << "Flow steady state: step "
            << simulation.steady_state_step()
            << " (t = "
            << simulation.steady_state_time()
            << " s)\n";
    }
    else
    {
        std::cout
            << "Flow steady state: not reached\n";
    }

    return 0;
}
