#include "counterflow/Combustion.hpp"

#include <algorithm>
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


ReactiveState advance_reaction_state(
    const ReactiveState& state,
    double density,
    double heat_capacity,
    double dt
)
{
    if (heat_capacity <= 0.0)
    {
        throw std::invalid_argument(
            "Heat capacity must be strictly positive."
        );
    }

    if (dt <= 0.0)
    {
        throw std::invalid_argument(
            "Chemical time step must be strictly positive."
        );
    }

    const ReactionRates rates =
        compute_methane_reaction_rates(
            state.ch4,
            state.o2,
            state.temperature,
            density
        );

    return {
        state.ch4
            + dt * rates.ch4 / density,

        state.o2
            + dt * rates.o2 / density,

        state.h2o
            + dt * rates.h2o / density,

        state.co2
            + dt * rates.co2 / density,

        state.temperature
            + dt
            * rates.heat_release
            / (density * heat_capacity)
    };
}


ChemicalSubcycling compute_chemical_subcycling(
    double hydro_dt,
    double maximum_temperature,
    double density
)
{
    if (hydro_dt <= 0.0)
    {
        throw std::invalid_argument(
            "Hydrodynamic time step must be strictly positive."
        );
    }

    if (maximum_temperature <= 0.0)
    {
        throw std::invalid_argument(
            "Maximum temperature must be strictly positive."
        );
    }

    if (density <= 0.0)
    {
        throw std::invalid_argument(
            "Density must be strictly positive."
        );
    }

    const double q_max =
        MethaneChemistry::pre_exponential_factor
        * (
            density * 0.5
            / MethaneChemistry::molar_mass_ch4
        )
        * std::pow(
            density * 0.1
            / MethaneChemistry::molar_mass_o2,
            2.0
        )
        * std::exp(
            -MethaneChemistry::activation_temperature
            / maximum_temperature
        );

    const double target_chemical_dt =
        1.0e3 / (10.0 * q_max);

    const double required_substeps =
        hydro_dt / target_chemical_dt;

    const std::size_t substeps =
        std::max<std::size_t>(
            1,
            static_cast<std::size_t>(
                std::ceil(required_substeps)
            )
        );

    return {
        substeps,
        hydro_dt / static_cast<double>(substeps)
    };
}

} // namespace counterflow
