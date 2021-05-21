# Mammals-restrict-alpine-plants

This repository contains the code used for analyses in the following publication:

Lynn, J.S., T.E.X. Miller, J.A. Rudgers. Accepted. Mammalian herbivores restrict the altitudinal range limits of alpine plants. Ecology Letters.

All the data used in the code can be found at the Environmental Data Initiative (EDI):

Lynn, J.S., J.A. Rudgers, and T.E. Miller. 2021. Mammalian herbivores restrict the altitudinal range limits of three alpine grass species (transplant and herbivore exclusion experiment and demographic data from natural populations), West Elk Mountains, Colorado, USA 2015-2018 ver 2. Environmental Data Initiative. https://doi.org/10.6073/pasta/193a9609b5ff5cec2690b3ac67b57c82

Once the data is downloaded from EDI, house it in the "data_lives_here" folder. 
Some smaller data files used in the analyses are also housed in the GitHub repository.

The pieces of code are structured as follows:

1) Analyses of herbivore damage in the experiment by focal species in the "Species_DamagMods_UngHerb.R" files.
2) Analyses of fitness metrics in the experiment by focal species in the "Species_FitnessMods_ungHerb.R" files.
3) Individual code for demographic analyses of each of the three focal species that 
  a) estimates vital rate functions
  b) apply experimental treatment effects to vital rate parameters
  c) estimates lambda with matrix projection models
  d) performs life table response experiment analyses
this code is is in the "Species_DemoMods_UngHerb.R" files.

Finally, code for the allometric equations used to predict focal species biomass is in the "AllometricEquationsFocalSppBiomass.R" file.