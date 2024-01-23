This is a ropository for the force field of the hydrogel hFF03 with chromophore and possible glycan additions
It mainly exists to act as a code repository for the paper: How chromophore labels shape the structure and dynamics of a peptide hydrogel (Currently in Prerelease)

In this repository you can finde:
1. A modified Amber99-SB-ILDN forcefield for GROMACS that contains the ortho and para aminobenzoic acid as a residue type as well as linker and basic GLYCAN paramters for a K17 modification (Currently Unpublished)
2. Julia Code for Analysing persistence length and self assembly for hFF03 and its variant. 
3. Julia Code for a spatially resolved self diffusion coefficient of water around the hFF03 coiled-coil
4. Julia Code for to create the plots used in the Paper. Contains matplotlib over pygui.
5. Example files for the starting configurations and necessary mdp and topology files used in the Paper. Also contains short example trajectories that work 

While the codes have some documentation it is reletivly sparse. I tried to name all variables and function in a sensible way and hopefully selfexplanatory


The main Code is in Julia but several codes use Python matplotlib loaded through julia.
The input is handled by the package Chemfiles and the codes can technically handle all file formats that Chemfiles can handle.
It uses .gro files for the topology. These are read and filtered to create mask lists. It might be problematic to use a different file format here.

Code 2A) Persistence Length: Has a parser and can be used from the command line. The parser is explained in the header of the File.
Code 2B) Selfassembly Analyse:  Change the name of the topology and trajectory file in the code and run. You can change the Number of Cores available. This Code is moderatly RAM intensive and might nor run on every PC.
         The Results are saved in jld2 format. This can be read by any julia code which loads the JLD2 package. There is an example vvaluation under EXAMPLEFiles
Code 3A) Diffusion around hff03:  Change the name of the topology and trajectory file in the code and run. You can change the batchsize aswell to see the influence but it is not advised. Necessary trajectory files are really huge.
     3B) Diffusion in Water: 

Code 4) Uses the JLD2 files obtaines by Code 2B

Example Files will be added when the manuscript is added. 


