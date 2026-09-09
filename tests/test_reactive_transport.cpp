#include <cmath>
#include <iostream>

#include "counterflow/Grid.hpp"
#include "counterflow/ReactiveTransport.hpp"

bool approximately_equal(
    double a,
    double b,
    double tolerance
)
{
    return std::abs(a - b) < tolerance;
}

int main()
{
    const counterflow::Grid2D grid(
        1.0,
        1.0,
        5,
        5
    );

    counterflow::Field2D ch4(
        grid.nx,
        grid.ny,
        0.5
    );

    counterflow::Field2D o2(
        grid.nx,
        grid.ny,
        0.1
    );

    counterflow::Field2D h2o(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D co2(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D temperature(
        grid.nx,
        grid.ny,
        1000.0
    );

    // Zero velocity isolates chemistry from advection.
    counterflow::Field2D u(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D v(
        grid.nx,
        grid.ny,
        0.0
    );

    const double density = 1.1614;
    const double heat_capacity = 1200.0;

    // Zero diffusivity isolates reaction as well.
    const double diffusivity = 0.0;

    const double hydro_dt = 1.0e-8;

    counterflow::ReactiveTransportStepper stepper(
        grid,
        density,
        heat_capacity,
        diffusivity,
        hydro_dt,
        1000.0
    );

    if (stepper.substeps() < 1)
    {
        std::cerr
            << "Reactive stepper has no chemical substeps.\n";

        return 1;
    }

    const double ch4_before =
        ch4(2, 2);

    const double o2_before =
        o2(2, 2);

    const double h2o_before =
        h2o(2, 2);

    const double co2_before =
        co2(2, 2);

    const double temperature_before =
        temperature(2, 2);

    stepper.advance(
        ch4,
        o2,
        h2o,
        co2,
        temperature,
        u,
        v,
        grid,
        true
    );

    if (!(ch4(2, 2) < ch4_before))
    {
        std::cerr
            << "CH4 was not consumed.\n";

        return 1;
    }

    if (!(o2(2, 2) < o2_before))
    {
        std::cerr
            << "O2 was not consumed.\n";

        return 1;
    }

    if (!(h2o(2, 2) > h2o_before))
    {
        std::cerr
            << "H2O was not produced.\n";

        return 1;
    }

    if (!(co2(2, 2) > co2_before))
    {
        std::cerr
            << "CO2 was not produced.\n";

        return 1;
    }

    if (!(temperature(2, 2) > temperature_before))
    {
        std::cerr
            << "Temperature did not increase.\n";

        return 1;
    }

    // Boundary conditions from the counterflow configuration.
    if (!approximately_equal(
            o2(0, 0),
            0.21,
            1.0e-12
        ))
    {
        std::cerr
            << "Bottom oxygen inlet is incorrect.\n";

        return 1;
    }

    if (!approximately_equal(
            ch4(0, grid.ny - 1),
            1.0,
            1.0e-12
        ))
    {
        std::cerr
            << "Top methane inlet is incorrect.\n";

        return 1;
    }

    return 0;
}
