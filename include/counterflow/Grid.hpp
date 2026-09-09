#pragma once

#include <cstddef>
#include <vector>

namespace counterflow
{

struct Grid2D
{
    double lx;
    double ly;

    std::size_t nx;
    std::size_t ny;

    double dx;
    double dy;

    std::vector<double> x;
    std::vector<double> y;

    Grid2D(
        double lx_,
        double ly_,
        std::size_t nx_,
        std::size_t ny_
    );

    [[nodiscard]]
    std::size_t index(std::size_t i, std::size_t j) const;
};


class Field2D
{
public:
    Field2D(
        std::size_t nx,
        std::size_t ny,
        double initial_value = 0.0
    );

    [[nodiscard]]
    std::size_t nx() const;

    [[nodiscard]]
    std::size_t ny() const;

    double& operator()(std::size_t i, std::size_t j);

    const double& operator()(std::size_t i, std::size_t j) const;

private:
    std::size_t nx_;
    std::size_t ny_;

    std::vector<double> data_;
};

} // namespace counterflow
