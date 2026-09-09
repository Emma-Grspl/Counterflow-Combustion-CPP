#include "counterflow/Grid.hpp"

#include <stdexcept>

namespace counterflow
{

Grid2D::Grid2D(
    double lx_,
    double ly_,
    std::size_t nx_,
    std::size_t ny_
)
    : lx(lx_),
      ly(ly_),
      nx(nx_),
      ny(ny_),
      dx(lx_ / static_cast<double>(nx_ - 1)),
      dy(ly_ / static_cast<double>(ny_ - 1)),
      x(nx_ * ny_),
      y(nx_ * ny_)
{
    if (nx < 2 || ny < 2)
    {
        throw std::invalid_argument(
            "Grid dimensions must contain at least two points."
        );
    }

    for (std::size_t j = 0; j < ny; ++j)
    {
        for (std::size_t i = 0; i < nx; ++i)
        {
            const std::size_t k = index(i, j);

            x[k] = static_cast<double>(i) * dx;
            y[k] = static_cast<double>(j) * dy;
        }
    }
}


std::size_t Grid2D::index(
    std::size_t i,
    std::size_t j
) const
{
    if (i >= nx || j >= ny)
    {
        throw std::out_of_range("Grid index out of range.");
    }

    return i + j * nx;
}


Field2D::Field2D(
    std::size_t nx,
    std::size_t ny,
    double initial_value
)
    : nx_(nx),
      ny_(ny),
      data_(nx * ny, initial_value)
{
    if (nx_ == 0 || ny_ == 0)
    {
        throw std::invalid_argument(
            "Field dimensions must be strictly positive."
        );
    }
}


std::size_t Field2D::nx() const
{
    return nx_;
}


std::size_t Field2D::ny() const
{
    return ny_;
}


double& Field2D::operator()(
    std::size_t i,
    std::size_t j
)
{
    if (i >= nx_ || j >= ny_)
    {
        throw std::out_of_range("Field index out of range.");
    }

    return data_[i + j * nx_];
}


const double& Field2D::operator()(
    std::size_t i,
    std::size_t j
) const
{
    if (i >= nx_ || j >= ny_)
    {
        throw std::out_of_range("Field index out of range.");
    }

    return data_[i + j * nx_];
}

} // namespace counterflow
