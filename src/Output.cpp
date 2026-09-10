#include "counterflow/Output.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace counterflow
{

void write_simulation_csv(
    const Simulation& simulation,
    const std::filesystem::path& output_path
)
{
    const std::filesystem::path parent =
        output_path.parent_path();

    if (!parent.empty())
    {
        std::filesystem::create_directories(parent);
    }

    std::ofstream file(output_path);

    if (!file)
    {
        throw std::runtime_error(
            "Failed to open output file: "
            + output_path.string()
        );
    }

    file << std::scientific
         << std::setprecision(17);

    file
        << "x,y,u,v,p,Y_N2,Y_CH4,Y_O2,Y_H2O,Y_CO2,T\n";

    const Grid2D& grid =
        simulation.grid();

    const Field2D& u =
        simulation.u();

    const Field2D& v =
        simulation.v();

    const Field2D& pressure =
        simulation.pressure();

    const Field2D& nitrogen =
        simulation.nitrogen();

    const Field2D& ch4 =
        simulation.ch4();

    const Field2D& o2 =
        simulation.o2();

    const Field2D& h2o =
        simulation.h2o();

    const Field2D& co2 =
        simulation.co2();

    const Field2D& temperature =
        simulation.temperature();

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            const std::size_t k =
                grid.index(i, j);

            file
                << grid.x[k] << ','
                << grid.y[k] << ','
                << u(i, j) << ','
                << v(i, j) << ','
                << pressure(i, j) << ','
                << nitrogen(i, j) << ','
                << ch4(i, j) << ','
                << o2(i, j) << ','
                << h2o(i, j) << ','
                << co2(i, j) << ','
                << temperature(i, j)
                << '\n';
        }
    }

    if (!file)
    {
        throw std::runtime_error(
            "Failed while writing output file: "
            + output_path.string()
        );
    }
}

} // namespace counterflow
