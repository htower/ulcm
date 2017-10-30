#include "core.h"

model::abstract::abstract(uint32_t n)
    : n(n)
    , fc(0)
    , gc(0)
{
}

method::abstract::abstract(const std::string& name)
    : epsilon(1e-4)
    , iter_max(1e7)
    , name(name)
    , iter(0)
    , time_i(0.0)
    , log_file(nullptr)
{
}

method::abstract::~abstract()
{
}

void
method::abstract::optimize(model::abstract& model, point& p)
{
    model.f (p);
    model.df(p);

    iter = 0;

    printf("%6s %7s %6s  %-23s %-23s %-8s %-8s\n",
           "method", "time", "iter", "f", "nrm2_g", "fc", "gc");
    print_status(model, p);
    printf("\n");

    log_file = fopen((name + ".log").c_str(), "w");
    fprintf(log_file, "%s;%s;%s;%s;%s;%s;%s\n", "method", "time", "iter", "f", "nrm2_g", "fc", "gc");
    log_status(model, p);

    time_0 = chrono::now();
    optimization_body(model, p);

    model.f (p);
    model.df(p);
    printf("\n");
    print_status(model, p);
    printf("\n\n");

    fclose(log_file);
}

bool
method::abstract::log_and_check(model::abstract& model, const point& p)
{
    ++iter;
    time_i = chrono::s(time_0);

    print_status(model, p);
    log_status(model, p);

    if (p.f <= epsilon * 5.0)
    {
        return true;
    }

    if (iter >= iter_max)
    {
        return true;
    }

    return false;
}

void
method::abstract::print_status(model::abstract& model, const point& p) const
{
    printf("%6s %7.3f %6u % 20.16e % 20.16e %8u %8u\r",
           name.c_str(),
           time_i,
           iter,
           p.f,
           p.gn,
           model.fc,
           model.gc);
}

void
method::abstract::log_status(model::abstract& model, const point& p) const
{
    fprintf(log_file, "%s;%e;%u;%20.16e;%20.16e;%u;%u\n",
           name.c_str(),
           time_i,
           iter,
           p.f,
           p.gn,
           model.fc,
           model.gc);
}

point::point(uint32_t n)
    : x(n)
    , g(n)
    , f(0.0)
    , gn(0.0)
{
}

point::point(const model::abstract& model)
    : point(model.n)
{
}

void
point::swap(point& other)
{
    x.swap(other.x);
    g.swap(other.g);

    std::swap(f,  other.f);
    std::swap(gn, other.gn);
}

ods::abstract::abstract(model::abstract& model)
    : model(model)
    , w(model)
{
}

ods::abstract::~abstract()
{
}

double
ods::abstract::f_alpha(const point& p, double alpha)
{
    w.x = p.x - alpha * p.g;
    model.f(w);
    return w.f;
}
