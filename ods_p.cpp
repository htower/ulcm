#include "ods_p.h"

namespace ods
{

parabolic_1p::parabolic_1p(model::abstract& model)
    : abstract(model)
{
}

parabolic_1p::~parabolic_1p()
{
}

double
parabolic_1p::search(const point& p)
{
    const double alpha_par = 1e0;
    const double df_0 = -blas::dot(p.g, p.g);

    double alpha_min = -0.5 * (df_0 * alpha_par * alpha_par) /
                              (f_alpha(p, alpha_par) - p.f - df_0 * alpha_par);

    return std::isnormal(alpha_min) ? alpha_min : 0.0;
}

}
