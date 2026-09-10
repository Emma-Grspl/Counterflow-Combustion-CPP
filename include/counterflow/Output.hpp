#pragma once

#include <filesystem>

#include "counterflow/Simulation.hpp"

namespace counterflow
{

void write_simulation_csv(
    const Simulation& simulation,
    const std::filesystem::path& output_path
);

} // namespace counterflow
