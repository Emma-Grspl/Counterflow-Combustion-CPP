#include "counterflow/ReactiveTransport.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "counterflow/Transport.hpp"

namespace counterflow
{

void initialize_reactive_fields(
    Field2D& ch4,
    Field2D& o2,
    Field2D& h2o,
    Field2D& co2,
    Field2D& temperature,
    const Grid2D& grid
)
{
    const bool dimensions_match =
        ch4.nx() == grid.nx &&
        ch4.ny() == grid.ny &&
        o2.nx() == grid.nx &&
        o2.ny() == grid.ny &&
        h2o.nx() == grid.nx &&
        h2o.ny() == grid.ny &&
        co2.nx() == grid.nx &&
        co2.ny() == grid.ny &&
        temperature.nx() == grid.nx &&
        temperature.ny() == grid.ny;

    if (!dimensions_match)
    {
        throw std::invalid_argument(
            "Reactive fields must match the grid."
        );
    }

    const std::size_t first_quarter =
        grid.nx / 4;

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            ch4(i, j) = 0.0;
            o2(i, j) = 0.0;
            h2o(i, j) = 0.0;
            co2(i, j) = 0.0;

            const std::size_t k =
                grid.index(i, j);

            const double y =
                grid.y[k];

            if (y >= 0.00075 &&
                y < 0.00125)
            {
                temperature(i, j) = 1000.0;
            }
            else
            {
                temperature(i, j) = 300.0;
            }
        }
    }

    // Oxygen inlet: bottom boundary.
    for (std::size_t i = 0;
         i < first_quarter;
         ++i)
    {
        o2(i, 0) = 0.21;
    }

    // Methane inlet: top boundary.
    for (std::size_t i = 0;
         i < first_quarter;
         ++i)
    {
        ch4(i, grid.ny - 1) = 1.0;
    }
}



void update_nitrogen_from_mass_closure(
    Field2D& nitrogen,
    const Field2D& ch4,
    const Field2D& o2,
    const Field2D& h2o,
    const Field2D& co2,
    const Grid2D& grid
)
{
    const bool dimensions_match =
        nitrogen.nx() == grid.nx &&
        nitrogen.ny() == grid.ny &&
        ch4.nx() == grid.nx &&
        ch4.ny() == grid.ny &&
        o2.nx() == grid.nx &&
        o2.ny() == grid.ny &&
        h2o.nx() == grid.nx &&
        h2o.ny() == grid.ny &&
        co2.nx() == grid.nx &&
        co2.ny() == grid.ny;

    if (!dimensions_match)
    {
        throw std::invalid_argument(
            "Species fields must match the grid."
        );
    }

    constexpr double tolerance = 1.0e-10;

    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        for (std::size_t i = 0; i < grid.nx; ++i)
        {
            const double reactive_sum =
                ch4(i, j)
                + o2(i, j)
                + h2o(i, j)
                + co2(i, j);

            if (reactive_sum < -tolerance ||
                reactive_sum > 1.0 + tolerance)
            {
                throw std::runtime_error(
                    "Reactive mass fractions violate "
                    "the mixture closure."
                );
            }

            nitrogen(i, j) =
                std::clamp(
                    1.0 - reactive_sum,
                    0.0,
                    1.0
                );
        }
    }
}


void apply_reactive_boundary_conditions(
    Field2D& ch4,
    Field2D& o2,
    Field2D& h2o,
    Field2D& co2,
    Field2D& temperature,
    const Grid2D& grid
)
{
    const bool dimensions_match =
        ch4.nx() == grid.nx &&
        ch4.ny() == grid.ny &&
        o2.nx() == grid.nx &&
        o2.ny() == grid.ny &&
        h2o.nx() == grid.nx &&
        h2o.ny() == grid.ny &&
        co2.nx() == grid.nx &&
        co2.ny() == grid.ny &&
        temperature.nx() == grid.nx &&
        temperature.ny() == grid.ny;

    if (!dimensions_match)
    {
        throw std::invalid_argument(
            "Reactive fields must match the grid."
        );
    }

    const std::size_t first_quarter =
        grid.nx / 4;

    const std::size_t half =
        grid.nx / 2;

    // Left and right boundaries:
    // zero normal gradient.
    for (std::size_t j = 0; j < grid.ny; ++j)
    {
        ch4(0, j) = ch4(1, j);
        o2(0, j) = o2(1, j);
        h2o(0, j) = h2o(1, j);
        co2(0, j) = co2(1, j);
        temperature(0, j) = temperature(1, j);

        ch4(grid.nx - 1, j) =
            ch4(grid.nx - 2, j);

        o2(grid.nx - 1, j) =
            o2(grid.nx - 2, j);

        h2o(grid.nx - 1, j) =
            h2o(grid.nx - 2, j);

        co2(grid.nx - 1, j) =
            co2(grid.nx - 2, j);

        temperature(grid.nx - 1, j) =
            temperature(grid.nx - 2, j);
    }

    // Bottom and top boundaries:
    // zero-gradient except at imposed inlets.
    for (std::size_t i = 0; i < grid.nx; ++i)
    {
        ch4(i, 0) = ch4(i, 1);
        o2(i, 0) = o2(i, 1);
        h2o(i, 0) = h2o(i, 1);
        co2(i, 0) = co2(i, 1);
        temperature(i, 0) =
            temperature(i, 1);

        ch4(i, grid.ny - 1) =
            ch4(i, grid.ny - 2);

        o2(i, grid.ny - 1) =
            o2(i, grid.ny - 2);

        h2o(i, grid.ny - 1) =
            h2o(i, grid.ny - 2);

        co2(i, grid.ny - 1) =
            co2(i, grid.ny - 2);

        temperature(i, grid.ny - 1) =
            temperature(i, grid.ny - 2);

        // Oxygen inlet at the bottom.
        if (i < first_quarter)
        {
            o2(i, 0) = 0.21;
        }

        // Methane inlet at the top.
        if (i < first_quarter)
        {
            ch4(i, grid.ny - 1) = 1.0;
        }

        // Cold inlet regions.
        if (i < half)
        {
            temperature(i, 0) = 300.0;
            temperature(i, grid.ny - 1) = 300.0;
        }
    }
}


ReactiveTransportStepper::ReactiveTransportStepper(
    const Grid2D& grid,
    double density,
    double heat_capacity,
    double diffusivity,
    double hydro_dt,
    double maximum_temperature
)
    : nx_(grid.nx),
      ny_(grid.ny),
      density_(density),
      heat_capacity_(heat_capacity),
      diffusivity_(diffusivity),
      hydro_dt_(hydro_dt),
      subcycling_(
          compute_chemical_subcycling(
              hydro_dt,
              maximum_temperature,
              density
          )
      ),
      ch4_buffer_(grid.nx, grid.ny, 0.0),
      o2_buffer_(grid.nx, grid.ny, 0.0),
      h2o_buffer_(grid.nx, grid.ny, 0.0),
      co2_buffer_(grid.nx, grid.ny, 0.0),
      temperature_buffer_(grid.nx, grid.ny, 0.0)
{
    if (density_ <= 0.0)
    {
        throw std::invalid_argument(
            "Density must be strictly positive."
        );
    }

    if (heat_capacity_ <= 0.0)
    {
        throw std::invalid_argument(
            "Heat capacity must be strictly positive."
        );
    }

    if (diffusivity_ < 0.0)
    {
        throw std::invalid_argument(
            "Diffusivity cannot be negative."
        );
    }
}


void ReactiveTransportStepper::advance(
    Field2D& ch4,
    Field2D& o2,
    Field2D& h2o,
    Field2D& co2,
    Field2D& temperature,
    const Field2D& u,
    const Field2D& v,
    const Grid2D& grid,
    bool enable_heat_release
)
{
    const bool dimensions_match =
        grid.nx == nx_ &&
        grid.ny == ny_ &&
        ch4.nx() == nx_ &&
        ch4.ny() == ny_ &&
        o2.nx() == nx_ &&
        o2.ny() == ny_ &&
        h2o.nx() == nx_ &&
        h2o.ny() == ny_ &&
        co2.nx() == nx_ &&
        co2.ny() == ny_ &&
        temperature.nx() == nx_ &&
        temperature.ny() == ny_ &&
        u.nx() == nx_ &&
        u.ny() == ny_ &&
        v.nx() == nx_ &&
        v.ny() == ny_;

    if (!dimensions_match)
    {
        throw std::invalid_argument(
            "Reactive transport fields must match the grid."
        );
    }

    double maximum_temperature =
        temperature(0, 0);

    for (std::size_t j = 0; j < ny_; ++j)
    {
        for (std::size_t i = 0; i < nx_; ++i)
        {
            maximum_temperature =
                std::max(
                    maximum_temperature,
                    temperature(i, j)
                );
        }
    }

    subcycling_ =
        compute_chemical_subcycling(
            hydro_dt_,
            maximum_temperature,
            density_
        );

    for (std::size_t substep = 0;
         substep < subcycling_.substeps;
         ++substep)
    {
        const double dt =
            subcycling_.chemical_dt;

        // ----------------------------------------------------
        // 1. Advection-diffusion
        // ----------------------------------------------------

        advance_species_transport(
            ch4,
            u,
            v,
            ch4_buffer_,
            grid,
            diffusivity_,
            dt
        );

        advance_species_transport(
            o2,
            u,
            v,
            o2_buffer_,
            grid,
            diffusivity_,
            dt
        );

        advance_species_transport(
            h2o,
            u,
            v,
            h2o_buffer_,
            grid,
            diffusivity_,
            dt
        );

        advance_species_transport(
            co2,
            u,
            v,
            co2_buffer_,
            grid,
            diffusivity_,
            dt
        );

        advance_scalar_transport(
            temperature,
            u,
            v,
            temperature_buffer_,
            grid,
            diffusivity_,
            dt
        );

        // Swap buffers instead of copying complete fields.
        std::swap(ch4, ch4_buffer_);
        std::swap(o2, o2_buffer_);
        std::swap(h2o, h2o_buffer_);
        std::swap(co2, co2_buffer_);
        std::swap(
            temperature,
            temperature_buffer_
        );

        // ----------------------------------------------------
        // 2. Local chemical reaction
        // ----------------------------------------------------

        for (std::size_t j = 1;
             j < grid.ny - 1;
             ++j)
        {
            for (std::size_t i = 1;
                 i < grid.nx - 1;
                 ++i)
            {
                const ReactiveState state{
                    ch4(i, j),
                    o2(i, j),
                    h2o(i, j),
                    co2(i, j),
                    temperature(i, j)
                };

                const ReactiveState reacted =
                    advance_reaction_state(
                        state,
                        density_,
                        heat_capacity_,
                        dt
                    );

                ch4(i, j) =
                    std::clamp(
                        reacted.ch4,
                        0.0,
                        1.0
                    );

                o2(i, j) =
                    std::clamp(
                        reacted.o2,
                        0.0,
                        1.0
                    );

                h2o(i, j) =
                    std::clamp(
                        reacted.h2o,
                        0.0,
                        1.0
                    );

                co2(i, j) =
                    std::clamp(
                        reacted.co2,
                        0.0,
                        1.0
                    );

                if (enable_heat_release)
                {
                    temperature(i, j) =
                        reacted.temperature;
                }
            }
        }

        // ----------------------------------------------------
        // 3. Boundary conditions
        // ----------------------------------------------------

        apply_reactive_boundary_conditions(
            ch4,
            o2,
            h2o,
            co2,
            temperature,
            grid
        );
    }
}


std::size_t ReactiveTransportStepper::substeps() const
{
    return subcycling_.substeps;
}


double ReactiveTransportStepper::chemical_dt() const
{
    return subcycling_.chemical_dt;
}

} // namespace counterflow
