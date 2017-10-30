#include "ufgm.h"

namespace method
{

ufgm::ufgm()
    : abstract("UFGM")
{
}

ufgm::~ufgm()
{
}

void
ufgm::optimization_body(model::abstract& model, point& p)
{
    dvector& x_0 = p.x;
    dvector  v_k(model.n);
    dvector  vvv(model.n);
    vvv.setZero();

    dvector dyx(model.n);

    point   y_k = p;
    point   x_kp1(model);
    point   y_kp1(model);
    dvector z_kp1(model.n);

    double alpha_k   = 0.0;
    double alpha_kp1 = alpha_k;
    double l_k       = 1.0;
    double l_kp1     = l_k;

    while (true)
    {
        l_kp1 = 0.5 * l_k;

        while (true)
        {
            v_k = x_0 - vvv;

            alpha_kp1 = 0.5 / l_kp1 + sqrt(0.25 / (l_kp1 * l_kp1) + alpha_k * alpha_k * l_k / l_kp1);
            double tau_k = 1.0 / (alpha_kp1 * l_kp1);

            x_kp1.x = tau_k * v_k + (1.0 - tau_k) * y_k.x;

            model.f (x_kp1);
            model.df(x_kp1);

            z_kp1 = v_k - alpha_kp1 * x_kp1.g;
            y_kp1.x = tau_k * z_kp1 + (1.0 - tau_k) * y_k.x;

            model.f(y_kp1);

            dyx = y_kp1.x - x_kp1.x;

            double cond = x_kp1.f - y_kp1.f;
            cond += blas::dot(x_kp1.g, dyx);
            cond += 0.5 * l_kp1 * blas::dot(dyx, dyx);
            cond += 0.5 * (tau_k * epsilon);

            if (cond >= 0.0)
            {
                break;
            }
            else
            {
                l_kp1 *= 2.0;
            }
        }

        l_k     = l_kp1;
        alpha_k = alpha_kp1;

        y_k.swap(y_kp1);
        vvv += alpha_kp1 * x_kp1.g;

        model.df(y_k);

        if (log_and_check(model, y_k))
        {
            break;
        }
    }

    p = y_k;
}

}
