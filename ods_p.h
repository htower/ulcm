#pragma once

#include "core.h"

namespace ods
{

struct parabolic_1p : public abstract
{
    parabolic_1p(model::abstract& model);
    ~parabolic_1p() override;

    double
    search(const point& p) override;
};

}
