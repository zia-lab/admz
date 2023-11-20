# Revisiting ADMZ

Repo of the old admz code with changes to make it run in modern machines.

Staring point from code (v 1.0) as downloaded from [here](https://elsevier.digitalcommonsdata.com/datasets/t6xsgx957b/1).

To compile or use must load `lapack` in CCV using `module load lapack`.

To facilitate generating the Vk matrix elements and then diagonalizing the matrix the bash script `runit.sh` should be copied to the folder where the data for an ion is located and run like `./runit.sh NUM_STATES NUM_ELECTRONS`.

Article that accompanies the code:
- Edvardsson, Sverker, and Daniel Åberg. “An Atomic Program for Energy Levels of Equivalent Electrons: Lanthanides and Actinides.” Computer Physics Communications 133, no. 2–3 (2001): 396–406.
