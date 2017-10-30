#pragma once

#include "core.h"


namespace method
{

struct ufgm : public abstract
{
    ufgm();
    ~ufgm() override;

protected:
    void
    optimization_body(model::abstract& model, point& p) override;
};

}
