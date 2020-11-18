# Randomized 1Q GST 1024 experiments

This directory contains files for randomizing 1Q experiments consisting
of sequences of Clifford group operations.

`randomize.jl` will take a sequence file (a CSV file, representing
Clifford group operations by a integer index between 1 and 24) and
generating a sequence corresponding to the insertion of (uniformly)
random Pauli gates. Strictly speaking these Pauli gates are not
independent, as we construct the sequences such that the measured
population is the same as the unrandomized population (at least for
ideal gates).

This randomization can proceed in two distinct ways: with compilation
and without compilation. Compilation refers to the combination of one
of the origonal Clifford group operations with the random Pauli group
operation, so that the experimental sequences length is unchanged. If
compilation is disabled, a sequence of length `n` is transformed into
a sequence of length `2n+1`.

`verify.jl` will compare a randomized sequences with an unrandomized
sequence, and verify that the sequences are equal up to a random
180deg rotation about Z.

