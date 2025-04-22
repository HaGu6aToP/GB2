#pragma once

#include "headers.h"
#include "types.h"

extern "C"{

    struct F4Result{
        Basis basis;
        ulong len;
    };

    typedef struct F4Result F4Result;

    F4Result F4(const Basis F, ulong npoly, const Field field, const PolynomRing ctx);

}