# MTBBuilder
A workflow of scripts to build a 13 protofilament microtubule from the component tubulin monomers. The workflow builds a topology psf and pdb using psfgen[] in VMD[]. There are also scripts to build microtubule bundles in different geometric bundles from a linear arrangement to hexagonal closest packed structure. Finally there is a script to convert the bundle to a format necessary to calculate scattering curves of the . The code is primarily written in TCL but contains python3 code for the final step. 

The structures is based on the 3JAL pdb structure. 
