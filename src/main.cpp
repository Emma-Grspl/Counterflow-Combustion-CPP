#include <iostream>

#include "counterflow/Config.hpp"

int main()
{
    const counterflow::SimulationConfig config;

    std::cout << "Counterflow Combustion C++ Solver\n";
    std::cout << "Grid: "
              << config.nx << " x " << config.ny << '\n';

    std::cout << "Time step: "
              << config.dt() << " s\n";

    return 0;
}
