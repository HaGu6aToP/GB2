# GB

This library implements the Buchberg algorithm, the improved Buchberger algorithm and F4 for ideals over finite fields $p < 2^{32}$. This code depends on FLINT, GLib and GBLA.

- FLINT is a C library for doing number theory, freely available under the GNU Lesser General Public License version 3 or later.
- GBLA is an open source (GPLv2) C library for linear algebra specialized for eliminating matrices generated during Gröbner basis computations in algorithms like F4 or F5. 
- GLib is a low-level library that extends the features provided by the standard C language libc library.

<!-- ## Dependence -->

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

|system   |F4   |Buchberger|macaulay2 Buchberger|macaulay2 F4|
|---------|-----|----------|--------------------|------------|
|cyclic4  |0.004|0.00006   |0.00013             |0.007       |
|cyclic5  |0.013|0.0025    |0.004               |0.06        |
|cyclic6  |0.04 |0.08      |0.06                |0.18        |
|cyclic7  |1.4  |16        |9                   |1.9         |
|cyclic8  |46   |>30m      |439                 |14          |
|katsura7 |0.04 |0.3       |0.8                 |0.15        |
|katsura8 |1.3  |2.8       |9                   |0.3         |
|katsura9 |8    |22        |101                 |1.1         |
|katsura10|81   |>30m      |>30m                |5           |
|katsura11|769  |>30m      |>30m                |32          |
|alea6    |4    |7         |349                 |1.7         |
|bayes148 |>30m |120       |12                  |88          |
|eco12    |207  |>30m      |>30m                |37          |
|noon9    |>30m |176       |>30m                |38          |
|reimer7  |285  |>30m      |>30m                |9           |
|schwarz11|9    |5         |64                  |1.1         |


$$
    p = 10^8 + 7
$$
|system   |F4   |Buchberger|macaulay2 Buchberger|macaulay2 F4|
|---------|-----|----------|--------------------|------------|
|cyclic4  |0.005|0.00006   |0.0004              |0.0003      |
|cyclic5  |0.016|0.00256   |0.004               |0.004       |
|cyclic6  |0.04 |0.08      |0.06                |0.08        |
|cyclic7  |1.5  |17        |9                   |12          |
|cyclic8  |48   |>30m      |283                 |458         |
|katsura7 |0.19 |0.4       |0.8                 |0.9         |
|katsura8 |1.3  |3         |9                   |9           |
|katsura9 |8    |23        |92                  |103         |
|katsura10|82   |>30m      |>30m                |1380        |
|katsura11|773  |>30m      |>30m                |>30m        |
|alea6    |4    |7         |37                  |49          |
|bayes148 |>30m |121       |13                  |14          |
|eco12    |213  |>30m      |>30m                |>30m        |
|noon9    |>30m |176       |>30m                |>10m        |
|reimer7  |290  |>30m      |>30m                |>10m        |
|schwarz11|10   |5         |10                  |13          |



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
