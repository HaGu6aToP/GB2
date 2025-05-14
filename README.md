# GB2

This library implements the Buchberg algorithm, the improved Buchberger algorithm, and F4 for ideals over finite fields $p < 2^{64}$. This code depends on FLINT, GLib and GBLA.

- FLINT is a C library for doing number theory, freely available under the GNU Lesser General Public License version 3 or later.
- GBLA is an open source (GPLv2) C library for linear algebra specialized for eliminating matrices generated during Gröbner basis computations in algorithms like F4 or F5. 
- GLib is a low-level library that extends the features provided by the standard C language libc library.

## Dependence

Algorithm F4 uses the functions of reducing the sparse matrix to row echelon form from GBLA. You can download it from <a href='https://hpac.imag.fr/gbla/'>link</a> (version 0.2 from 2016). Then you need to go to the mapping.c file and change lines 758 

``` M->pos[i][M->rwidth[i]] = map->pc[i];``` 

and 837 

```M->pos[i][M->rwidth[i]] = map->npc_rev[map->npiv+j];``` 

to 

```M->pos[i][M->rwidth[i]]   = i;``` 

and 

```M->pos[i][M->rwidth[i]]   = map->npiv+j;``` 

(Or I can provide a compiled version by me if someone sees this and asks for it :D). Then follow the instructions in GBLA.

