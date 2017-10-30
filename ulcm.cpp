#include "ulcm.h"

namespace method
{

ulcm::ulcm()
    : abstract("ULCM")
{
}

ulcm::~ulcm()
{
}

void
ulcm::optimization_body(model::abstract& model, point& p)
{
    dvector dz(model.n);

    point   y_k = p;
    dvector z_k = p.x;
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
            alpha_kp1 = 0.5 / l_kp1 + sqrt(0.25 / (l_kp1 * l_kp1) + alpha_k * alpha_k * l_k / l_kp1);
            double tau_k = 1.0 / (alpha_kp1 * l_kp1);

            x_kp1.x = tau_k * z_k + (1.0 - tau_k) * y_k.x;

            model.f (x_kp1);
            model.df(x_kp1);
            double h_kp1 = model.ods().search(x_kp1);

            y_kp1.x = x_kp1.x - h_kp1 * x_kp1.g;
            model.f(y_kp1);

            z_kp1 = z_k - alpha_kp1 * x_kp1.g;

            dz = z_k - z_kp1;

            double cond = alpha_kp1 * blas::dot(x_kp1.g, dz);
            cond -= 0.5 * blas::dot(dz, dz);
            cond -= alpha_kp1 * alpha_kp1 * l_kp1 * (x_kp1.f - y_kp1.f + 0.5 * tau_k * epsilon);

            if (cond <= 0.0)
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
        z_k.swap(z_kp1);

        model.df(y_k);

        if (log_and_check(model, y_k))
        {
            break;
        }
    }

    p = y_k;
}

}
