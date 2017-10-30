#pragma once

#include "core.h"

namespace method
{

struct ulcm : public method::abstract
{
    ulcm();
    ~ulcm() override;

protected:
    void
    optimization_body(model::abstract& model, point & p) override;
};

}
