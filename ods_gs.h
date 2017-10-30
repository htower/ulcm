#pragma once

#include "core.h"

namespace ods
{

struct golden_section : public abstract
{
    golden_section(model::abstract& model);
    ~golden_section() override;

    double
    search(const point& p) override;

    double alpha_0;
    double tol;
};

}
