#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>

#include "counterflow/Grid.hpp"
#include "counterflow/PressurePoisson.hpp"


int main()
{
    constexpr std::size_t nx = 50;
    constexpr std::size_t ny = 50;

    constexpr double lx = 0.002;
    constexpr double ly = 0.002;

    constexpr std::size_t warmup_solves = 100;
    constexpr std::size_t solves_per_batch = 1000;
    constexpr std::size_t measured_batches = 7;

    const counterflow::Grid2D grid(
        lx,
        ly,
        nx,
        ny
    );

    counterflow::Field2D rhs(
        nx,
        ny,
        0.0
    );

    counterflow::Field2D pressure(
        nx,
        ny,
        0.0
    );

    // --------------------------------------------------------
    // Deterministic non-zero RHS.
    //
    // Boundary values remain zero because those rows encode
    // the pressure boundary conditions.
    // --------------------------------------------------------

    for (
        std::size_t j = 1;
        j < ny - 1;
        ++j
    )
    {
        for (
            std::size_t i = 1;
            i < nx - 1;
            ++i
        )
        {
            const double x =
                static_cast<double>(i)
                / static_cast<double>(nx - 1);

            const double y =
                static_cast<double>(j)
                / static_cast<double>(ny - 1);

            rhs(i, j) =
                std::sin(3.141592653589793 * x)
                * std::cos(2.0 * 3.141592653589793 * y);
        }
    }

    using Clock =
        std::chrono::steady_clock;

    // --------------------------------------------------------
    // Measure construction + sparse factorization separately.
    // --------------------------------------------------------

    const auto construction_start =
        Clock::now();

    counterflow::PressurePoissonSolver solver(
        grid
    );

    const auto construction_end =
        Clock::now();

    const double construction_seconds =
        std::chrono::duration<double>(
            construction_end
            - construction_start
        ).count();

    std::cout
        << std::fixed
        << std::setprecision(9);

    std::cout
        << "C++ sparse pressure benchmark\n"
        << "=============================\n"
        << "Grid                 : "
        << nx
        << " x "
        << ny
        << '\n'
        << "Unknowns             : "
        << solver.system_size()
        << '\n'
        << "Matrix nonzeros      : "
        << solver.nonzero_count()
        << '\n'
        << "Construction/factor. : "
        << construction_seconds
        << " s\n\n";

    // --------------------------------------------------------
    // Warm-up.
    // --------------------------------------------------------

    for (
        std::size_t n = 0;
        n < warmup_solves;
        ++n
    )
    {
        solver.solve(
            rhs,
            pressure
        );
    }

    std::cout
        << "Warm-up solves       : "
        << warmup_solves
        << '\n'
        << "Solves per batch     : "
        << solves_per_batch
        << '\n'
        << "Measured batches     : "
        << measured_batches
        << "\n\n";

    std::vector<double> per_solve_times;

    per_solve_times.reserve(
        measured_batches
    );

    double checksum = 0.0;

    for (
        std::size_t batch = 0;
        batch < measured_batches;
        ++batch
    )
    {
        const auto start =
            Clock::now();

        for (
            std::size_t n = 0;
            n < solves_per_batch;
            ++n
        )
        {
            solver.solve(
                rhs,
                pressure
            );
        }

        const auto end =
            Clock::now();

        const double batch_seconds =
            std::chrono::duration<double>(
                end - start
            ).count();

        const double per_solve =
            batch_seconds
            / static_cast<double>(
                solves_per_batch
            );

        per_solve_times.push_back(
            per_solve
        );

        // Force observable use of the result.
        checksum +=
            pressure(nx / 2, ny / 2);

        std::cout
            << "batch "
            << batch + 1
            << ": "
            << batch_seconds
            << " s total, "
            << per_solve * 1.0e6
            << " us/solve\n";
    }

    // --------------------------------------------------------
    // Basic statistics.
    // --------------------------------------------------------

    double sum = 0.0;

    for (
        const double value :
        per_solve_times
    )
    {
        sum += value;
    }

    const double mean =
        sum
        / static_cast<double>(
            per_solve_times.size()
        );

    std::vector<double> sorted =
        per_solve_times;

    std::sort(
        sorted.begin(),
        sorted.end()
    );

    const double median =
        sorted[
            sorted.size() / 2
        ];

    double variance = 0.0;

    for (
        const double value :
        per_solve_times
    )
    {
        const double delta =
            value - mean;

        variance +=
            delta * delta;
    }

    variance /=
        static_cast<double>(
            per_solve_times.size() - 1
        );

    const double stddev =
        std::sqrt(variance);

    const double estimated_2999 =
        median * 2999.0;

    std::cout
        << "\nSummary\n"
        << "-------\n"
        << "mean solve          : "
        << mean * 1.0e6
        << " us\n"
        << "median solve        : "
        << median * 1.0e6
        << " us\n"
        << "std                 : "
        << stddev * 1.0e6
        << " us\n"
        << "2999 solve estimate : "
        << estimated_2999
        << " s\n"
        << "checksum            : "
        << checksum
        << '\n';

    return 0;
}
