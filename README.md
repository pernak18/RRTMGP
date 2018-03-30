# First Iteration

We are rewriting the RRTMGP cloud code, and I figured it wouldn't hurt to do this separate designing and coding in a different version control repo. The original, untouched cloud code with Robert Pincus notes are in mo_cloud_optics_pade.F90, and my revisions will be written in pernak_mo_cloud_optics_pade.F90.

# Second Iteration

All of the previous code has been relocated to a subdirectory called iter1.

Starting over, kinda. After speaking with Robert (on 29-Mar-2018), we decided that the eventual goal is to expand the cloud optics class from by-band to by-g-point using random sampling (to map cloud optical properties to g-point spectral dimension) and overlap (to structure random numbers in the vertical direction, forcing similar layers to be close to each other) similar to what's happening in $RRTMGP/extensions/cloud_sampling.F90. Analogous to gas optics, initialization will happen by receiving all inputs arguments (gas optics and level/layer T and P state variables) and storing them in the data (tables), all of which will be returned via an object of either a 1-scaler, 2-stream, or n-stream class.

The test driver for cloud optical properties still needs to generated, as do the data (which can just be Garand 1 with cloud parameters like effective radius, liquid and ice water paths, and cloud fraction). have a look at the validation scripts and Frank's add_cloud_optics (?) to see how this can be done.

Some relevant clarifications from Robert:

 - the relevant classes for the proposed optical properties class are:
   1) *spectral discretization ($RRTMGP/rte/mo_spectral_disc.F90)*: gas optics
   2) *optical properties ($RRTMGP/rte/mo_optical_props.F90)*: radiative transfer dimensions
- For 2), every spectral calculation is done independently

- the optical properties class is extended by the array class, which are 3-D structures in which the data are stored
- no data are in the two base classes -- the data are provided by the user and stored at initialization via the array class
- as the objects "grow", new data (e.g., tau) needs to be added to the existing data
- but optical property and array objects are abstract in that they cannot be instantiated; rather, the array objects are extended by "concrete" 1-scalar, 2-stream, and n-stream subclasses that can be instantiated once an analytical solution (e.g., n-stream) is derived from the input arguments. these subclasses are what will be returned

