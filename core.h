#pragma once

#include "Eigen/Dense"

#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <chrono>
#include <cmath>
#include <vector>

using dvector = Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor>;

namespace blas
{

inline double
dot(const dvector& v1, const dvector& v2)
{
    return v1.dot(v2);
}

inline double
nrm2(const dvector& v)
{
    return v.stableNorm();
}

}

namespace chrono
{

using point = std::chrono::time_point<std::chrono::system_clock>;

inline point
now()
{
    return std::chrono::system_clock::now();
}

inline double
ms(const point& start, const point& end)
{
    return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start).count();
}

inline double
ms(const point& start)
{
    return ms(start, now());
}

inline double
s(const point& start, const point& end)
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

inline double
s(const point& start)
{
    return s(start, now());
}

}

struct point;

namespace ods
{

struct abstract;

}

namespace model
{

struct abstract
{
    abstract(uint32_t n);

    virtual void
    f(point& p) = 0;

    virtual void
    df(point& p) = 0;

    virtual ods::abstract&
    ods() = 0;

    uint32_t n;

    uint32_t fc;
    uint32_t gc;
};

}

namespace method
{

struct abstract
{
    abstract(const std::string& name);
    virtual ~abstract();

    void
    optimize(model::abstract& model, point& p);

    double epsilon;

    uint32_t iter_max;

protected:
    virtual void
    optimization_body(model::abstract& model, point& p) = 0;

    bool
    log_and_check(model::abstract& model, const point& p);

    void
    print_status(model::abstract& model, const point& p) const;

    void
    log_status(model::abstract& model, const point& p) const;

    std::string name;
    uint32_t iter;

    chrono::point time_0;
    double time_i;

    FILE* log_file;
};

}

struct point
{
    point(uint32_t n);

    point(const model::abstract& model);

    void
    swap(point& other);

    dvector x;
    dvector g;

    double f;
    double gn;
};

namespace ods
{

struct abstract
{
    abstract(model::abstract& model);
    virtual ~abstract();

    virtual double
    search(const point& p) = 0;

protected:
    double
    f_alpha(const point& p, double alpha);

    model::abstract& model;
    point w;
};

}
