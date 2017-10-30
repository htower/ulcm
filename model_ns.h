#pragma once

#include "core.h"

namespace model
{

struct non_smooth : public abstract
{
    non_smooth(uint32_t n);

    void
    f(point& p) override;

    void
    df(point& p) override;

    ods::abstract&
    ods() override;

    double mu;
};

}
