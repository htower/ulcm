#pragma once

#include "core.h"

namespace model
{

struct smooth : public abstract
{
    smooth(uint32_t n);

    void
    f(point& p) override;

    void
    df(point& p) override;

    ods::abstract&
    ods() override;
};

}
