#pragma once

#include "headers.h"
#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

    struct F4Result{
        Basis basis;
        ulong len;
    };

    typedef struct F4Result F4Result;

    F4Result F4(const Basis F, ulong npoly, const Field field, const PolynomRing ctx);

#ifdef __cplusplus
}
#endif