#include "model_ns.h"

#include "ods_gs.h"

namespace model
{

non_smooth::non_smooth(uint32_t n)
    : abstract(n)
    , mu(1e-1)
{
}

void
non_smooth::f(point& p)
{
    double x_max = p.x[0];

    p.f = 0.0;
    for (uint32_t i = 0; i < n; ++i)
    {
        p.f += p.x[i] * p.x[i];
        x_max = std::max(x_max, p.x[i]);
    }

    p.f *= mu;
    p.f += x_max;

    ++fc;
}

void
non_smooth::df(point& p)
{
    uint32_t i_max = 0;
    double   x_max = p.x[i_max];

    for (uint32_t i = 0; i < n; ++i)
    {
        p.g[i] = 2.0 * mu * p.x[i];

        if (p.x[i] > x_max)
        {
            x_max = p.x[i];
            i_max = i;
        }
    }

    p.g[i_max] += 1.0;
    p.gn = blas::nrm2(p.g);

    ++gc;
}

ods::abstract&
non_smooth::ods()
{
    static ods::golden_section o(*this);

    return o;
}

}
