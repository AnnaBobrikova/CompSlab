The Figs. corresponding to Figs 1-2 in Salmi et al. 2021 can be produced as follows:

Fig 1. Left:
python3 beaming_plot.py Comp_e100mu9_ ang figs/beaming_angX.pdf


Fig 1. Right:
python3 poldeg_plot.py Comp_e100mu9_ ang figs/poldeg_angX.pdf


Fig 2. Left:
python3 beaming_plot.py Comp_e100mu9_ ene figs/beaming_specX.pdf


Fig 2. Right:
python3 poldeg_plot.py Comp_e100mu9_ ene figs/poldeg_specX.pdf

To create Figs corresponding to 3-4 in Salmi et al. 2021:
Run driver.py in order to save pulses for 1 or 2 hot spots (separately or in combination if the spots are antipodal). 
A couple of example pulses are already saved in the 'pulses' directory.
Then run plot_stokes_model.py (and manually check that name of the input files are correct there and switch between antipodal/non-antipodal flags if necessary)



