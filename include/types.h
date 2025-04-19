#pragma once
#include "headers.h"

typedef fq_nmod_mpoly_struct* Polynom;
typedef fq_nmod_mpoly_struct** Basis;
typedef fq_nmod_mpoly_ctx_struct* PolynomRing;
typedef fq_nmod_ctx_struct* Field;

struct Pair{
    ulong first;
    ulong second;
};

struct SPair{
    Polynom poly;
    ulong first;
    ulong second;
};

struct Buchberger_result{
    Basis basis;
    ulong len;
};

struct thread_buchberger_data_t{
    ulong* finished_threads;
    Polynom S_poly;
    Basis F;
    ulong npoli;
    PolynomRing ctx;
    int completed;
};

struct thread_buchberger_data_v2_t{
    int executed;
    int running;
    int pause;
    ulong pair;
    ulong* executed_threads;
    GArray* F; 
    Polynom S_poly;
    PolynomRing ctx;
    struct Pair pPair;
};

struct F4Pair{
    Polynom lcm;
    Polynom t_f;
    Polynom f;
    Polynom t_g;
    Polynom g;
};

struct F4PairProjection{
    Polynom t;
    Polynom f;
};

typedef struct Pair Pair;
typedef struct SPair SPair;
typedef struct Buchberger_result Buchberger_result;
typedef struct thread_buchberger_data_t thread_buchberger_data_t;
typedef struct thread_buchberger_data_v2_t thread_buchberger_data_v2_t;
typedef struct F4Pair F4Pair;
typedef struct F4PairProjection F4PairProjection;