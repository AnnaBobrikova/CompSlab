# Atmosphere-thomsonization
Computing intensity and polarization degree of emergent radiation from electron-scattering dominated atmosphere. 

Code structure
---------------------------------------

Model for Compton slab:
* cs.jl -  Most updated model for obtaining Stokes I and Q from a Compton slab

Most updated pulse profile modelling code:
* ixpe_tools/polpulse.py -  Main model for polarized pulse-spectral computation (pulse profile modelling)
* ixpe_tools/driver.py - Simple tool to run the previous function

Previous code versions:
* cs.py -  Main model for polarized pulse-spectral computation
* cs_ts2_func.py - Almost same as previous but as a function used e.g. for mcmc fitting. 
* driver_paeqs.py - Simple tool to run the previous function
* emcee_cs_ts.py - Fitting PA profiles using emcee. 
... more to be added ...
