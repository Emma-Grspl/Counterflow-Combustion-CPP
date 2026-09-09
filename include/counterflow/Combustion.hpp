#pragma once

namespace counterflow
{

struct ReactionRates
{
    double ch4;
    double o2;
    double h2o;
    double co2;
    double heat_release;
};


struct MethaneChemistry
{
    // Arrhenius parameters
    static constexpr double activation_temperature = 10000.0;
    static constexpr double pre_exponential_factor = 1.1e8;

    // Molar masses [kg/mol]
    static constexpr double molar_mass_ch4 = 16.04e-3;
    static constexpr double molar_mass_o2  = 31.99e-3;
    static constexpr double molar_mass_h2o = 18.01e-3;
    static constexpr double molar_mass_co2 = 44.00e-3;

    // Stoichiometric coefficients
    static constexpr double nu_ch4 = -1.0;
    static constexpr double nu_o2  = -2.0;
    static constexpr double nu_h2o =  2.0;
    static constexpr double nu_co2 =  1.0;

    // Formation enthalpies [J/mol]
    static constexpr double h_ch4 = -74.9e3;
    static constexpr double h_o2  = 0.0;
    static constexpr double h_h2o = -241.818e3;
    static constexpr double h_co2 = -393.52e3;
};


[[nodiscard]]
ReactionRates compute_methane_reaction_rates(
    double y_ch4,
    double y_o2,
    double temperature,
    double density
);

} // namespace counterflow
