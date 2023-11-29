# Revisiting ADMZ

Repo of the old admz code with changes to make it run in modern machines.

Staring point from code (v 1.0) as downloaded from [here](https://elsevier.digitalcommonsdata.com/datasets/t6xsgx957b/1).

To compile or use must load `lapack` in CCV using `module load lapack`.

Best way of running this is using the `lanthanide.py` script which on its own takes care of a few things that before had to be done manually. These convenience include compiling the binary, and creating the input files that it requires.

Article that accompanies the code:
- Edvardsson, Sverker, and Daniel Åberg. “An Atomic Program for Energy Levels of Equivalent Electrons: Lanthanides and Actinides.” Computer Physics Communications 133, no. 2–3 (2001): 396–406.
