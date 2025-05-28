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

    // Алгоритм F4
    F4Result F4(const Basis F, ulong npoly, const Field field, const PolynomRing ctx);

    // Редуцирование матрицу к верхне треугольному виду
    void F4_poly_reduce(GArray* F_ref, const GArray* F, const GArray* F_monoms, const Field field, const PolynomRing);


#ifdef __cplusplus
}
#endif
