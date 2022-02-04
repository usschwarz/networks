# ReaDDy 2 scripts for "A particle-based computational model to analyse remodelling of the red blood cell cytoskeleton during malaria infections "

The following scripts have been written by Julia Jäger for ReaDDy 2 to simulate
remodelling of the actin-spectrin cytoskeleton by the malaria parasite.

Several python packages are needed and ReaDDy needs to be installed as 
described here: https://readdy.github.io/installation.html

## Script 1: shear runs of the cytoskeleton

For shearing a network containing actin filaments, a series of Readdy simulations 
was set up and run with the following script: 

```bash
python simulation-shear-actin-network.py
```

## Script 2: analysis of actin lengths based on simulation runs

Once several networks were simulated with different concentrations, one can run 
the following script to extract the lengths of all the contained actin filaments:

```bash
python analysis-of-actin-dynamics.py
```

## Script 3: simulation of KAHRP

The simulations containing KAHRP were set up according to the following script:

```bash
python simulation-shear-actin-network.py
```
