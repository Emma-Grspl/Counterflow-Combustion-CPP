#include <cmath>
#include <iostream>

#include "counterflow/Combustion.hpp"

bool approximately_equal(double a, double b, double tolerance)
{
    return std::abs(a - b) < tolerance;
}

int main()
{
    const double density = 1.1614;

    const counterflow::ReactionRates rates =
        counterflow::compute_methane_reaction_rates(
            0.5,
            0.1,
            1000.0,
            density
        );

    // Reactants must be consumed.
    if (!(rates.ch4 < 0.0))
    {
        std::cerr << "Methane should be consumed.\n";
        return 1;
    }

    if (!(rates.o2 < 0.0))
    {
        std::cerr << "Oxygen should be consumed.\n";
        return 1;
    }

    // Products must be formed.
    if (!(rates.h2o > 0.0))
    {
        std::cerr << "Water should be produced.\n";
        return 1;
    }

    if (!(rates.co2 > 0.0))
    {
        std::cerr << "Carbon dioxide should be produced.\n";
        return 1;
    }

    // Combustion must release heat.
    if (!(rates.heat_release > 0.0))
    {
        std::cerr << "Methane combustion should release heat.\n";
        return 1;
    }

    // Stoichiometric molar-rate relations.
    const double methane_molar_rate =
        rates.ch4
        / counterflow::MethaneChemistry::molar_mass_ch4;

    const double oxygen_molar_rate =
        rates.o2
        / counterflow::MethaneChemistry::molar_mass_o2;

    if (!approximately_equal(
            oxygen_molar_rate,
            2.0 * methane_molar_rate,
            1.0e-8
        ))
    {
        std::cerr
            << "Incorrect CH4/O2 stoichiometric relation.\n";

        return 1;
    }


    // ========================================================
    // Test 2: local chemical time step
    // ========================================================

    const counterflow::ReactiveState state{
        0.5,    // CH4
        0.1,    // O2
        0.0,    // H2O
        0.0,    // CO2
        1000.0  // Temperature [K]
    };

    const double heat_capacity = 1200.0;
    const double chemical_dt = 1.0e-8;

    const counterflow::ReactiveState new_state =
        counterflow::advance_reaction_state(
            state,
            density,
            heat_capacity,
            chemical_dt
        );

    // Reactants decrease.
    if (!(new_state.ch4 < state.ch4))
    {
        std::cerr
            << "CH4 did not decrease during reaction.\n";

        return 1;
    }

    if (!(new_state.o2 < state.o2))
    {
        std::cerr
            << "O2 did not decrease during reaction.\n";

        return 1;
    }

    // Products increase.
    if (!(new_state.h2o > state.h2o))
    {
        std::cerr
            << "H2O did not increase during reaction.\n";

        return 1;
    }

    if (!(new_state.co2 > state.co2))
    {
        std::cerr
            << "CO2 did not increase during reaction.\n";

        return 1;
    }

    // Exothermic reaction raises temperature.
    if (!(new_state.temperature > state.temperature))
    {
        std::cerr
            << "Temperature did not increase during reaction.\n";

        return 1;
    }

    // The original state must remain unchanged.
    if (!approximately_equal(
            state.ch4,
            0.5,
            1.0e-12
        ))
    {
        std::cerr
            << "Input reactive state was modified.\n";

        return 1;
    }


    // ========================================================
    // Test 3: chemical subcycling
    // ========================================================

    const double hydro_dt = 1.0e-4;

    const counterflow::ChemicalSubcycling subcycling =
        counterflow::compute_chemical_subcycling(
            hydro_dt,
            1500.0,
            density
        );

    if (subcycling.substeps < 1)
    {
        std::cerr
            << "Chemical substep count must be positive.\n";

        return 1;
    }

    if (!(subcycling.chemical_dt > 0.0))
    {
        std::cerr
            << "Chemical time step must be positive.\n";

        return 1;
    }

    const double reconstructed_hydro_dt =
        static_cast<double>(subcycling.substeps)
        * subcycling.chemical_dt;

    if (!approximately_equal(
            reconstructed_hydro_dt,
            hydro_dt,
            1.0e-14
        ))
    {
        std::cerr
            << "Chemical substeps do not exactly cover "
               "the hydrodynamic time step.\n";

        return 1;
    }

    if (subcycling.chemical_dt > hydro_dt)
    {
        std::cerr
            << "Chemical time step cannot exceed "
               "hydrodynamic time step.\n";

        return 1;
    }

    return 0;
}
