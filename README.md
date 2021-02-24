# Flat Histogram Adsorption Simulation Analysis

## Interpret results from flat histogram adsorption simulations. 

This code is based on the NVT+W method discussed in many previous academic works. It forms the basis of calculations performed in my recent paper **(DOI: 10.1021/acs.jpcc.0c11082)**. 

Basically, the flat histogram Monte Carlo simulations for adsorption yield the probabilities of a system going from having N molecules to N+1 or N-1 molecules. 

From these, the favorability of each such loading marostate is determined. Then, histogram reweighting is used to calculate this distribution at required temperature and chemical potential conditions and adsorption behavior can be calculated. 