#include <algorithm>
#include <cmath>
#include <iostream>

#include "counterflow/Grid.hpp"
#include "counterflow/NavierStokes.hpp"

bool approximately_equal(double a, double b, double tolerance)
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

    // ========================================================
    // Test 1: predictor step
    // ========================================================

    counterflow::Field2D u_old(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D v_old(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D u_star(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D v_star(
        grid.nx,
        grid.ny,
        0.0
    );

    // u(x,y) = x
    // v(x,y) = 0

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            const std::size_t k =
                grid.index(i, j);

            u_old(i, j) =
                grid.x[k];
        }
    }

    const double dt = 0.1;
    const double nu = 0.01;

    counterflow::compute_intermediate_velocity(
        u_old,
        v_old,
        u_star,
        v_star,
        grid,
        dt,
        nu
    );

    // At the center:
    //
    // u = 0.5
    // du/dx = 1
    // Laplacian(u) = 0
    //
    // u* = 0.5 - 0.1 * 0.5 = 0.45

    if (!approximately_equal(
            u_star(2, 2),
            0.45,
            1.0e-12
        ))
    {
        std::cerr
            << "Unexpected u_star value: "
            << u_star(2, 2)
            << '\n';

        return 1;
    }

    if (!approximately_equal(
            v_star(2, 2),
            0.0,
            1.0e-12
        ))
    {
        std::cerr
            << "Unexpected v_star value.\n";

        return 1;
    }

    // ========================================================
    // Test 2: pressure correction
    // ========================================================

    counterflow::Field2D pressure(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D corrected_u(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D corrected_v(
        grid.nx,
        grid.ny,
        0.0
    );

    // Manufactured pressure:
    //
    // p(x,y) = 2x + 3y
    //
    // therefore:
    //
    // dp/dx = 2
    // dp/dy = 3

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            const std::size_t k =
                grid.index(i, j);

            pressure(i, j) =
                2.0 * grid.x[k]
                + 3.0 * grid.y[k];

            // Use simple constant intermediate velocities
            // for this independent test.
            u_star(i, j) = 10.0;
            v_star(i, j) = 20.0;
        }
    }

    const double rho = 2.0;
    const double correction_dt = 0.1;

    counterflow::correct_velocity(
        u_star,
        v_star,
        pressure,
        corrected_u,
        corrected_v,
        grid,
        rho,
        correction_dt
    );

    // dt / rho = 0.05
    //
    // u_new = 10 - 0.05 * 2 = 9.9
    // v_new = 20 - 0.05 * 3 = 19.85

    if (!approximately_equal(
            corrected_u(2, 2),
            9.9,
            1.0e-12
        ))
    {
        std::cerr
            << "Unexpected corrected u velocity: "
            << corrected_u(2, 2)
            << '\n';

        return 1;
    }

    if (!approximately_equal(
            corrected_v(2, 2),
            19.85,
            1.0e-12
        ))
    {
        std::cerr
            << "Unexpected corrected v velocity: "
            << corrected_v(2, 2)
            << '\n';

        return 1;
    }

    // Boundary values must remain equal to u_star/v_star.

    if (!approximately_equal(
            corrected_u(0, 2),
            10.0,
            1.0e-12
        ))
    {
        std::cerr
            << "Boundary u velocity was modified.\n";

        return 1;
    }


    // ========================================================
    // Test 3: complete fractional step
    // ========================================================

    counterflow::Field2D step_u_old(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D step_v_old(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D step_u_new(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D step_v_new(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D step_pressure(
        grid.nx,
        grid.ny,
        0.0
    );

    // Uniform velocity field.
    //
    // All spatial derivatives vanish, therefore:
    //
    // predictor        -> unchanged
    // pressure RHS     -> zero
    // pressure         -> zero
    // correction       -> unchanged

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            step_u_old(i, j) = 1.5;
            step_v_old(i, j) = -0.5;
        }
    }

    counterflow::NavierStokesStepper stepper(
        grid,
        1.0,
        0.01,
        0.1
    );

    stepper.advance(
        step_u_old,
        step_v_old,
        step_u_new,
        step_v_new,
        step_pressure,
        grid
    );

    if (!approximately_equal(
            step_u_new(2, 2),
            1.5,
            1.0e-12
        ))
    {
        std::cerr
            << "Complete step modified uniform u velocity.\n";

        return 1;
    }

    if (!approximately_equal(
            step_v_new(2, 2),
            -0.5,
            1.0e-12
        ))
    {
        std::cerr
            << "Complete step modified uniform v velocity.\n";

        return 1;
    }


    // ========================================================
    // Test 4: discrete projection must remove divergence
    // ========================================================

    counterflow::Field2D projection_u_star(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D projection_v_star(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D projection_rhs(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D projection_pressure(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D projection_u(
        grid.nx,
        grid.ny,
        0.0
    );

    counterflow::Field2D projection_v(
        grid.nx,
        grid.ny,
        0.0
    );

    // Deliberately divergent velocity field.
    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            const std::size_t k =
                grid.index(i, j);

            const double x =
                grid.x[k];

            const double y =
                grid.y[k];

            projection_u_star(i, j) =
                x * x + 0.4 * y;

            projection_v_star(i, j) =
                -0.3 * x + y * y;
        }
    }

    const double projection_rho = 1.3;
    const double projection_dt = 0.05;

    counterflow::compute_pressure_rhs(
        projection_u_star,
        projection_v_star,
        projection_rhs,
        grid,
        projection_rho,
        projection_dt
    );

    counterflow::PressurePoissonSolver
        projection_solver(grid);

    projection_solver.solve(
        projection_rhs,
        projection_pressure
    );

    counterflow::correct_velocity(
        projection_u_star,
        projection_v_star,
        projection_pressure,
        projection_u,
        projection_v,
        grid,
        projection_rho,
        projection_dt
    );

    double maximum_divergence = 0.0;

    for (std::size_t j = 1;
         j < grid.ny - 1;
         ++j)
    {
        for (std::size_t i = 1;
             i < grid.nx - 1;
             ++i)
        {
            const double du_dx =
                (projection_u(i, j)
                 - projection_u(i - 1, j))
                / grid.dx;

            const double dv_dy =
                (projection_v(i, j)
                 - projection_v(i, j - 1))
                / grid.dy;

            maximum_divergence =
                std::max(
                    maximum_divergence,
                    std::abs(du_dx + dv_dy)
                );
        }
    }

    if (maximum_divergence > 1.0e-10)
    {
        std::cerr
            << "Pressure projection left excessive "
               "discrete divergence: "
            << maximum_divergence
            << '\n';

        return 1;
    }

    return 0;
}
