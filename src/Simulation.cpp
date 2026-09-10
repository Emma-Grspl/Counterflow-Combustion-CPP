#include "counterflow/Simulation.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "counterflow/BoundaryConditions.hpp"

namespace counterflow
{

Simulation::Simulation(
    const SimulationConfig& config
)
    : config_(config),
      grid_(
          config.lx,
          config.ly,
          config.nx,
          config.ny
      ),
      u_(
          config.nx,
          config.ny,
          0.0
      ),
      v_(
          config.nx,
          config.ny,
          0.0
      ),
      u_next_(
          config.nx,
          config.ny,
          0.0
      ),
      v_next_(
          config.nx,
          config.ny,
          0.0
      ),
      pressure_(
          config.nx,
          config.ny,
          0.0
      ),
      nitrogen_(
          config.nx,
          config.ny,
          0.0
      ),
      ch4_(
          config.nx,
          config.ny,
          0.0
      ),
      o2_(
          config.nx,
          config.ny,
          0.0
      ),
      h2o_(
          config.nx,
          config.ny,
          0.0
      ),
      co2_(
          config.nx,
          config.ny,
          0.0
      ),
      temperature_(
          config.nx,
          config.ny,
          300.0
      ),
      navier_stokes_(
          grid_,
          config.rho,
          config.nu,
          config.dt()
      ),
      reactive_transport_(
          grid_,
          config.rho,
          config.cp,
          config.diffusivity,
          config.dt(),
          1000.0
      )
{
    initialize_vertical_velocity(
        v_,
        grid_
    );

    initialize_reactive_fields(
        ch4_,
        o2_,
        h2o_,
        co2_,
        temperature_,
        grid_
    );

    update_nitrogen_from_mass_closure(
        nitrogen_,
        ch4_,
        o2_,
        h2o_,
        co2_,
        grid_
    );
}


void Simulation::step()
{
    if (step_count_ >= config_.nt - 1)
    {
        throw std::runtime_error(
            "Simulation has already reached its final time."
        );
    }

    // ========================================================
    // 1. Hydrodynamics
    // ========================================================

    navier_stokes_.advance(
        u_,
        v_,
        u_next_,
        v_next_,
        pressure_,
        grid_
    );

    std::swap(
        u_,
        u_next_
    );

    std::swap(
        v_,
        v_next_
    );

    // ========================================================
    // 2. Reactive transport
    // ========================================================

    reactive_transport_.advance(
        ch4_,
        o2_,
        h2o_,
        co2_,
        temperature_,
        u_,
        v_,
        grid_,
        true
    );

    update_nitrogen_from_mass_closure(
        nitrogen_,
        ch4_,
        o2_,
        h2o_,
        co2_,
        grid_
    );

    ++step_count_;

    time_ =
        static_cast<double>(step_count_)
        * config_.dt();
}


void Simulation::run(
    std::size_t number_of_steps
)
{
    if (step_count_ + number_of_steps
        > config_.nt - 1)
    {
        throw std::invalid_argument(
            "Requested number of steps exceeds "
            "the configured final simulation time."
        );
    }

    for (std::size_t n = 0;
         n < number_of_steps;
         ++n)
    {
        step();
    }
}


std::size_t Simulation::step_count() const
{
    return step_count_;
}


double Simulation::time() const
{
    return time_;
}


double Simulation::max_temperature() const
{
    double maximum =
        temperature_(0, 0);

    for (std::size_t j = 0;
         j < grid_.ny;
         ++j)
    {
        for (std::size_t i = 0;
             i < grid_.nx;
             ++i)
        {
            maximum =
                std::max(
                    maximum,
                    temperature_(i, j)
                );
        }
    }

    return maximum;
}


const Grid2D& Simulation::grid() const
{
    return grid_;
}


const Field2D& Simulation::u() const
{
    return u_;
}


const Field2D& Simulation::v() const
{
    return v_;
}


const Field2D& Simulation::pressure() const
{
    return pressure_;
}


const Field2D& Simulation::nitrogen() const
{
    return nitrogen_;
}


const Field2D& Simulation::ch4() const
{
    return ch4_;
}


const Field2D& Simulation::o2() const
{
    return o2_;
}


const Field2D& Simulation::h2o() const
{
    return h2o_;
}


const Field2D& Simulation::co2() const
{
    return co2_;
}


const Field2D& Simulation::temperature() const
{
    return temperature_;
}

} // namespace counterflow
