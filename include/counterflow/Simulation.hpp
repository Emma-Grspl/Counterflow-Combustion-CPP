#pragma once

#include <cstddef>

#include "counterflow/Config.hpp"
#include "counterflow/Grid.hpp"
#include "counterflow/NavierStokes.hpp"
#include "counterflow/ReactiveTransport.hpp"

namespace counterflow
{

class Simulation
{
public:
    explicit Simulation(
        const SimulationConfig& config
    );

    void step();

    void run(
        std::size_t number_of_steps
    );

    [[nodiscard]]
    std::size_t step_count() const;

    [[nodiscard]]
    double time() const;

    [[nodiscard]]
    double max_temperature() const;

    [[nodiscard]]
    const Grid2D& grid() const;

    [[nodiscard]]
    const Field2D& u() const;

    [[nodiscard]]
    const Field2D& v() const;

    [[nodiscard]]
    const Field2D& pressure() const;

    [[nodiscard]]
    const Field2D& nitrogen() const;

    [[nodiscard]]
    const Field2D& ch4() const;

    [[nodiscard]]
    const Field2D& o2() const;

    [[nodiscard]]
    const Field2D& h2o() const;

    [[nodiscard]]
    const Field2D& co2() const;

    [[nodiscard]]
    const Field2D& temperature() const;

private:
    SimulationConfig config_;
    Grid2D grid_;

    Field2D u_;
    Field2D v_;

    Field2D u_next_;
    Field2D v_next_;

    Field2D pressure_;

    Field2D nitrogen_;
    Field2D nitrogen_next_;

    Field2D ch4_;
    Field2D o2_;
    Field2D h2o_;
    Field2D co2_;
    Field2D temperature_;

    NavierStokesStepper navier_stokes_;
    ReactiveTransportStepper reactive_transport_;

    std::size_t step_count_ = 0;
    double time_ = 0.0;
};

} // namespace counterflow
