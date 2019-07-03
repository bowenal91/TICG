# TICG_Gridless

This repository contains several versions of the Theoretically-Informed Coarse-Grained Model for polymer simulations. The standard version is located in src_grid. All other versions are modifications of my own to explore interesting systems.

src_grid: Standard TICG as described in Detcheverry et al (2008) in Macromolecules.
src: Gridless version with pair potentials. Incorporates rigid body and nematic interactions.
src_WLC: Incorporates a wormlike chain model by using Monte Carlo moves that enforce rigidity.
src_CNC; Incorporates both a three body potential and rigid cylinders to describe polymers grafted onto cellulose nanocrystals
src_TI_Distance and src_TI_Rotation: Same as the CNC model with the ability to calculate free energy as a function of distance and rotation, respectively using thermodynamic integration.

