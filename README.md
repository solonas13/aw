AW: Avoided Words
===

GNU GPLv3 License; Copyright (C) 2016 Jia Gao and Solon P. Pissis

To compile AW, please follow the instructions given in file INSTALL.

```
 Usage: aw <options>
 Standard (Mandatory):
  -a, --alphabet            <str>     `DNA' for nucleotide  sequences or `PROT'
                                      for protein  sequences. 
  -i, --input-file          <str>     (Multi)FASTA input filename.
  -o, --output-file         <str>     Output filename.
  -t, --threshold           <dbl>     The threshold for AWs.
 Optional:
  -k, --length              <int>     Fixed length for AWs (default: no fixed).
  -A, --absent              <int>     `1' to check also for absent avoided words
                                      or `0' otherwise (default: 0).
                                      This option cannot be used with `-c 1'.
  -c, --common              <int>     `1' to check for common words instead of
                                      avoided or `0' otherwise (default: 0).
                                      This option can only be used with `-k <int>'.
  -r, --reverse             <int>     `1' to check for the reverse complement or
                                      `0' otherwise (default: 0).
```

When publishing work that is based on the results from AW please cite:

```
Y. Almirantis, P. Charalampopoulos, J. Gao, C. S. Iliopoulos, M. Mohamed, S. P. Pissis, D. Polychronopoulos, "Optimal Computation of Avoided Words", in Algorithms in Bioinformatics: 16th International Workshop, WABI 2016, Aarhus, Denmark, August 22-24, 2016. Proceedings, M. Frith, N. C. Storm Pedersen, Eds., Springer International Publishing, 2016, pp. 1-13.
```
