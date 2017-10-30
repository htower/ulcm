#include "ods_gs.h"

namespace ods
{

golden_section::golden_section(model::abstract& model)
    : abstract(model)
    , alpha_0(1e-3)
    , tol(1e-3)
{
}

golden_section::~golden_section()
{
}

inline void
shift(double& a, double& b, double& c)
{
    a = b;
    b = c;
}

inline void
shift(double& a, double& b, double& c, double d)
{
    a = b;
    b = c;
    c = d;
}

double
golden_section::search(const point& p)
{
    static const double K  = 0.618034;
    static const double K1 = 1.618034;

    double a  = 0.0;
    double fa = p.f;
    double b  = alpha_0; // 1.0 / p.gn if gradient exists
    double fb = f_alpha(p, b);
    double c;
    double fc;

    if (fb > fa)
    {
        std::swap(a,  b );
        std::swap(fa, fb);
    }

    c  = b + K1 * (b - a);
    fc = f_alpha(p, c);

    while (fb > fc)
    {
        double u  = c + K1 * (c - b);
        double fu = f_alpha(p, u);

        shift(a,  b,  c,  u );
        shift(fa, fb, fc, fu);
    }

    double I   = K * (c - a);
    double a1  = b;
    double fa1 = fb;
    double c1  = a + I;
    double fc1 = f_alpha(p, c1);

    while (fabs(I) > tol)
    {
        I *= K;

        if (fa1 >= fc1)
        {
            shift(a,  a1,  c1 );
            shift(fa, fa1, fc1);

            c1  = a + I;
            fc1 = f_alpha(p, c1);
        }
        else
        {
            shift(c,  c1,  a1 );
            shift(fc, fc1, fa1);

            a1  = c - I;
            fa1 = f_alpha(p, a1);
        }
    }

    return (fa1 < fc1) ? a1 : c1;
}

}
