#include "model_s.h"

#include "ods_p.h"

namespace model
{

smooth::smooth(uint32_t n)
    : abstract(n)
{
}

void
smooth::f(point& p)
{
    p.f = 0.0;

    for (uint32_t i = 0; i < n; ++i)
    {
        p.f += (i + 1) * p.x[i] * p.x[i];
    }

    ++fc;
}

void
smooth::df(point& p)
{
    p.gn = 0.0;

    for (uint32_t i = 0; i < n; ++i)
    {
        p.g[i] = 2.0 * (i + 1) * p.x[i];
        p.gn += p.g[i] * p.g[i];
    }

    p.gn = sqrt(p.gn);
    ++gc;
}

ods::abstract&
smooth::ods()
{
    static ods::parabolic_1p o(*this);

    return o;
}

}
