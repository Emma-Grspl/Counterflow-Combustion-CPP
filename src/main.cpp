#include <charconv>
#include <cstddef>
#include <iostream>
#include <string_view>

#include "counterflow/Config.hpp"
#include "counterflow/Simulation.hpp"

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
    const counterflow::SimulationConfig config;

    std::size_t number_of_steps = 10;

    if (argc > 2)
    {
        std::cerr
            << "Usage: "
            << argv[0]
            << " [number_of_steps|all]\n";

        return 1;
    }

    if (argc == 2)
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
        << "\n\n";

    counterflow::Simulation simulation(
        config
    );

    simulation.run(
        number_of_steps
    );

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

    return 0;
}
