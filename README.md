# First Iteration

We are rewriting the RRTMGP cloud code, and I figured it wouldn't hurt to do this separate designing and coding in a different version control repo. The original, untouched cloud code with Robert Pincus notes are in mo_cloud_optics_pade.F90, and my revisions will be written in pernak_mo_cloud_optics_pade.F90.

# Second Iteration

## Notes From Robert

All of the previous code has been relocated to a subdirectory called iter1.

Starting over, kinda. After speaking with Robert (on 29-Mar-2018), we decided that the eventual goal is to expand the cloud optics class from by-band to by-g-point using random sampling (to map cloud optical properties to g-point spectral dimension) and overlap (to structure random numbers in the vertical direction, forcing similar layers to be close to each other) similar to what's happening in $RRTMGP/extensions/cloud_sampling.F90. Analogous to gas optics, initialization will happen by receiving all inputs arguments (gas optics and level/layer T and P state variables) and storing them in the data (tables), all of which will be returned via an object of either a 1-scalar, 2-stream, or n-stream class.

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

I've copied over all of the relevant RRTMGP code from the SVN repository for reference when designing the cloud optics class. We will first do this in Python, then port to Fortran 2003.

## Notes From Mike

 - cloud physical properties: water path, effective radius, cloud fraction
 - cloud optical properties: single scatter albedo (ssa or omega), asymmetry (g), optical depth (tau, extinction)
 - physical properties are provided by the user (e.g., GCM modelers), and we compute the corresponding optical properties are for all of our bands, then at our g-points. for now, the capability for computation onto both spectral domains exists just in case only one of the conversions (instead of both) is needed
 - we have two ways of doing this: with lookup table (LUT) "data" and with the Pade approximation
 - netCDF files for each are in the data/ subdirectory of Mike's cloud code development branch (https://aervcs-ext.aer.com/aer/rd/RRTMGP/branches/mji_dev-branch/data/)
 - the roughness dimension for ice comes from the Yang paper that Eli forwarded to Mike and me many months ago
 - column dimension in the cloud netCDF input/output files (for LUT only?) corresponds to the number of clouds considered in the computations, which is done in conjunction with a single clear sky profile
 - the "size regime" dimension exists because Mike had to fit the Pade curves to the LUT data, but this fitting did not work well when considering all possible effective radii. instead, Mike broke up the r_eff into multiple regimes (via trial and error) and then performed the fits for both spectral domains and water phases. each pade_sizereg_* array then contains 4 elements that represent the boundaries of each r_eff domain used in the fits. THE ARRAYS SHOULD NOT BE INPUT BY THE USER AND SHOULD BE KEPT IN THEIR CURRENT HARD-CODED STATE
 - Pade Coefficients and LUT values are also hard-coded, effectively. they will be provided in netCDF files that will be read by our code, and the user will not have the option of providing their own coefficient arrays
 - coefficient loading code is in $RRTMGP/test/util/src/io (read_cldpp  and read_cldop are just physical property and optical property readers, respectively)
 - class() and type() functions are "interchangeable" for our purposes
 - in the phase function array, asymmetry is the first moment
 - cloud particle sizes are between 2.5 and 21 microns; particles that are any bigger are considered precipitation

