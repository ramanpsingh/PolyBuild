# PolyBuild

Build custom Gromacs coordinate and topology files (.gro and .itp) for Martini/CG model of Pluronic surfactants.

Pluronics have a general structure of PEOx-PPOy-PEOx, where x is the number of EO monomers and y is the number of PO monomers. 

Usage:
python3 PolyBuild.py -peo [Number of PEO beads in one block] -ppo [Number of PPO beads in the central block] -bl [Bond length] -o [Output file name]
All fields are optional. .If an option is not used, a default values will be used.
The parameters are based on:
1. Perez-Sanchez, G., et al. (2019). Rationalizing the Phase Behavior of Triblock Copolymers through Experiments and Molecular Simulations. The Journal of Physical Chemistry C 123(34): 21224-21236. https://doi.org/10.1021/acs.jpcc.9b04099
2. Pérez-Sánchez, G., et al. (2021). Using coarse-grained molecular dynamics to understand the effect of ionic liquids on the aggregation of Pluronic copolymer solutions. Physical Chemistry Chemical Physics 23(10): 5824-5833. https://doi.org/10.1039/D0CP06572B
3. Hezaveh, S., et al. (2012). Understanding the Interaction of Block Copolymers with DMPC Lipid Bilayer Using Coarse-Grained Molecular Dynamics Simulations. The Journal of Physical Chemistry B 116(49): 14333-14345. https://doi.org/10.1021/jp306565e

For example:
python3 PolyBuild.py -peo 100 -ppo 50
This will generate a Martini CG model corresponding to PEO100-PPO50-PEO100. A default bondlength of 0.280 will be used. 

The script can be easily adapted for other linear polymers as well.

