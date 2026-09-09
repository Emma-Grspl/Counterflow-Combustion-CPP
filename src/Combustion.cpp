#include "counterflow/Combustion.hpp"

#include <cmath>
#include <stdexcept>

namespace counterflow
{

ReactionRates compute_methane_reaction_rates(
    double y_ch4,
    double y_o2,
    double temperature,
    double density
)
{
    if (temperature <= 0.0)
    {
        throw std::invalid_argument(
            "Temperature must be strictly positive."
        );
    }

    if (density <= 0.0)
    {
        throw std::invalid_argument(
            "Density must be strictly positive."
        );
    }

    const double methane_concentration =
        density * y_ch4
        / MethaneChemistry::molar_mass_ch4;

    const double oxygen_concentration =
        density * y_o2
        / MethaneChemistry::molar_mass_o2;

    const double arrhenius_factor =
        std::exp(
            -MethaneChemistry::activation_temperature
            / temperature
        );

    const double reaction_rate =
        MethaneChemistry::pre_exponential_factor
        * methane_concentration
        * oxygen_concentration
        * oxygen_concentration
        * arrhenius_factor;

    const double w_ch4 =
        MethaneChemistry::molar_mass_ch4
        * MethaneChemistry::nu_ch4
        * reaction_rate;

    const double w_o2 =
        MethaneChemistry::molar_mass_o2
        * MethaneChemistry::nu_o2
        * reaction_rate;

    const double w_h2o =
        MethaneChemistry::molar_mass_h2o
        * MethaneChemistry::nu_h2o
        * reaction_rate;

    const double w_co2 =
        MethaneChemistry::molar_mass_co2
        * MethaneChemistry::nu_co2
        * reaction_rate;

    const double heat_release =
        -(
            MethaneChemistry::h_ch4
            / MethaneChemistry::molar_mass_ch4
            * w_ch4

            + MethaneChemistry::h_o2
            / MethaneChemistry::molar_mass_o2
            * w_o2

            + MethaneChemistry::h_h2o
            / MethaneChemistry::molar_mass_h2o
            * w_h2o

            + MethaneChemistry::h_co2
            / MethaneChemistry::molar_mass_co2
            * w_co2
        );

    return {
        w_ch4,
        w_o2,
        w_h2o,
        w_co2,
        heat_release
    };
}

} // namespace counterflow
