#include "model_s.h"
#include "model_ns.h"
#include "ulcm.h"
#include "ufgm.h"

int
main(int, char**)
{
    const uint32_t n       = 1e3;
    const double   epsilon = 1e-4;

    model::smooth     model(n);
//     model::non_smooth model(n);

    dvector x0(n);
    x0.setConstant(10.0);

    method::ufgm ufgm;
    method::ulcm ulcm;
    point        p0(n);

    p0.x = x0;
    ufgm.epsilon = epsilon;
    ufgm.optimize(model, p0);

    p0.x = x0;
    ulcm.epsilon = epsilon;
    ulcm.optimize(model, p0);

    return 0;
}
