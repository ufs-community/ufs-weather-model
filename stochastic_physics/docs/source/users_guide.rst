Users Guide    
==================================================
The stochastic physics currently only works with the UFS-atmosphere model

Currently, 3 stochastic schemes are used operationally at NCEP/EMC: Stochastic Kinetic Energy Backscatter (SKEB; Berner et al., 2009), Stochastically Perturbed Physics Tendencies (SPPT; Palmer et al., 2009), and Specific Humidity perturbations (SHUM), which is inspired by Tompkins and Berner, 2008. In addition there is the ability to perturb certain land model/surface parameters (Gehne et al, 2019), and a cellular automata scheme (Bengtsson et al. 2019) which interacts directly with the convective parameterization.

SKEB adds wind perturbations to model state.  Perturbations are random in space/time, but amplitude is determined by a smoothed dissipation estimate provided by the dynamical core. 
Addresses errors in the dynamics  - more active in the mid-latitudes

SPPT multiplies the physics tendencies by a random number O [0,2] before updating the model state.  Addresses error in the physics parameterizations (either missing physics or unresolved subgrid processes). It is most active in boundary layer and convective regions

SHUM multiply the low-level specific humidity by a small random number each time-step. It attempts to address missing physics (cold pools, gust fronts), most active in convective regions

Land surface perturbations allow for land surface parameters such as Albedo, Soil Hydraulic Conductivity, LAI, and roughness lengths to vary in space. Addresses error in the land model and land-atmosphere interactions.  

Due to the model’s numerics, any stochastic perturbation needs to be correlated in space and time in order to have the desired effect of upscale growth of the perturbations. This is achieved by creating a random pattern that has a specified decorrelation length-scale and is a first order auto-regressive process AR(1) in time with a specified decorrelation time-scale.  (The CA random pattern generator also satisfies this condition)

Currently the Land surface perturbations and cellular automata are not supported at the workflow level.  

