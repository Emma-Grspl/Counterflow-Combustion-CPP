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


    // ========================================================
    // Test 1: physical initialization
    // ========================================================

    counterflow::initialize_reactive_fields(
        ch4,
        o2,
        h2o,
        co2,
        temperature,
        grid
    );

    if (!approximately_equal(
            o2(0, 0),
            0.21,
            1.0e-12
        ))
    {
        std::cerr
            << "Incorrect initial oxygen inlet.\n";

        return 1;
    }

    if (!approximately_equal(
            ch4(0, grid.ny - 1),
            1.0,
            1.0e-12
        ))
    {
        std::cerr
            << "Incorrect initial methane inlet.\n";

        return 1;
    }

    // Reset to a homogeneous reacting state for the
    // reactive-step test below.
    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            ch4(i, j) = 0.5;
            o2(i, j) = 0.1;
            h2o(i, j) = 0.0;
            co2(i, j) = 0.0;
            temperature(i, j) = 1000.0;
        }
    }

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

    const std::size_t initial_substeps =
        stepper.substeps();

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

    if (stepper.substeps() < initial_substeps)
    {
        std::cerr
            << "Adaptive chemical subcycling unexpectedly "
               "reduced the number of substeps.\n";

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


    // ========================================================
    // Complete inlet compositions
    // ========================================================

    const std::size_t bottom_air_i = 0;
    const std::size_t slow_n2_i = grid.nx / 4;

    // Bottom fast jet: air.
    if (!approximately_equal(
            o2(bottom_air_i, 0),
            0.21,
            1.0e-12
        ) ||
        !approximately_equal(
            ch4(bottom_air_i, 0),
            0.0,
            1.0e-12
        ) ||
        !approximately_equal(
            h2o(bottom_air_i, 0),
            0.0,
            1.0e-12
        ) ||
        !approximately_equal(
            co2(bottom_air_i, 0),
            0.0,
            1.0e-12
        ))
    {
        std::cerr
            << "Bottom air inlet composition is incorrect.\n";

        return 1;
    }

    // Top fast jet: pure methane.
    if (!approximately_equal(
            ch4(0, grid.ny - 1),
            1.0,
            1.0e-12
        ) ||
        !approximately_equal(
            o2(0, grid.ny - 1),
            0.0,
            1.0e-12
        ) ||
        !approximately_equal(
            h2o(0, grid.ny - 1),
            0.0,
            1.0e-12
        ) ||
        !approximately_equal(
            co2(0, grid.ny - 1),
            0.0,
            1.0e-12
        ))
    {
        std::cerr
            << "Top methane inlet composition is incorrect.\n";

        return 1;
    }

    // Slow bottom and top jets: pure nitrogen.
    // Reactive species must therefore all be zero.
    if (!approximately_equal(
            ch4(slow_n2_i, 0),
            0.0,
            1.0e-12
        ) ||
        !approximately_equal(
            o2(slow_n2_i, 0),
            0.0,
            1.0e-12
        ) ||
        !approximately_equal(
            h2o(slow_n2_i, 0),
            0.0,
            1.0e-12
        ) ||
        !approximately_equal(
            co2(slow_n2_i, 0),
            0.0,
            1.0e-12
        ) ||
        !approximately_equal(
            ch4(slow_n2_i, grid.ny - 1),
            0.0,
            1.0e-12
        ) ||
        !approximately_equal(
            o2(slow_n2_i, grid.ny - 1),
            0.0,
            1.0e-12
        ))
    {
        std::cerr
            << "Slow nitrogen inlet composition is incorrect.\n";

        return 1;
    }

    return 0;
}
