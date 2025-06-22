# GB2

This library implements the Buchberg algorithm, the improved Buchberger algorithm and F4 for ideals over finite fields $p < 2^{32}$. This code depends on FLINT, GLib and GBLA.

- FLINT is a C library for doing number theory, freely available under the GNU Lesser General Public License version 3 or later.
- GBLA is an open source (GPLv2) C library for linear algebra specialized for eliminating matrices generated during Gröbner basis computations in algorithms like F4 or F5. 
- GLib is a low-level library that extends the features provided by the standard C language libc library.

## Dependence

<!-- Algorithm F4 uses the functions of reducing the sparse matrix to row echelon form from GBLA. You can download it from <a href='https://hpac.imag.fr/gbla/'>link</a> (version 0.2 from 2016). Then you need to go to the **mapping.c** file and change (function reconstruct_matrix_block_no_multiline) lines 758 

``` M->pos[i][M->rwidth[i]] = map->pc[i];``` 

and 837 

```M->pos[i][M->rwidth[i]] = map->npc_rev[map->npiv+j];``` 

to 

```M->pos[i][M->rwidth[i]]   = i;``` 

and 

```M->pos[i][M->rwidth[i]]   = map->npiv+j;``` 

Also all "if" block in "write B part" on

```
M->rows[block_row_idx + min_range_blocks - j - 1][M->rwidth[block_row_idx + min_range_blocks - j - 1]]  = B->blocks[l][i].val[k+line_idx];

M->pos[block_row_idx + min_range_blocks - j - 1][M->rwidth[block_row_idx + min_range_blocks - j - 1]]   = map->npiv + k + start_idx;//map->npc_rev[map->npiv+k+start_idx];

M->rwidth[block_row_idx + min_range_blocks - j - 1]++;
```

And in "write D part"

```map->npc_rev[map->npiv+j]```

to

```map->npiv+j```

(Or I can provide a compiled version by me if someone sees this and asks for it :D). Then follow the installation instructions in GBLA. -->

## BenchMark

$$
    p = 30011
$$

|system   |F4   |buchberger|macaulay2|
|---------|-----|----------|---------|
|cyclic4  |0.004|0.00006   |0.00013  |
|cyclic5  |0.013|0.0025    |0.0035   |
|cyclic6  |0.04 |0.08      |0.06     |
|cyclic7  |1.4  |16        |9        |
|cyclic8  |46   |>10m      |439      |
|katsura7 |0.04 |0.3       |0.8      |
|katsura8 |1.3  |2.8       |9        |
|katsura9 |8    |22        |101      |
|katsura10|81   |-         |>30m     |
|alea6    |4    |7         |349      |
|bayes148 |>30m |120       |75       |
|eco12    |207  |-         |>30m     |
|noon9    |>30m |176       |-        |
|reimer7  |285  |-         |-        |
|schwarz11|9    |5         |64       |


<!-- |system |F4 (sec)|
|-------|--------|
|bayes148| >7200|
|cyclic7|0.055|
|cyclic8|0.37|
|cyclic9|2.068|
|eco12|0.140|
|jason210|memory error(>900)|
|katrusra7|0.036|
|katrusra8|0.055|
|katrusra9|0.024|
|katrusra10|0.158|
|katrusra11|0.312|
|katrusra12|0.561|
|mayr42|8.435|
|noon9|0.813|
|reimer7|0.478|
|reimer8|9.269|
|schwar11|0.502|
|yang1|>1200| -->
