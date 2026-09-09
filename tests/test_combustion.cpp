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

    return 0;
}
