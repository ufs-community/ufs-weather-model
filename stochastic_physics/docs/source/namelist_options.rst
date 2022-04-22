Stochastic Physics Namelist 
===========================

General options 
"""""""""""""""

.. csv-table::
   :header: "Option", "Description"
   :widths: 30, 50

   "NTRUNC", "Optional, Spectral resolution (e.g. T126) of random patterns, default is for model to determine proper truncation"
   "LAT_S", "Optional, number of latitude points for the gaussian grid  (must be even), default is for model to determine gaussian grid"
   "LON_S", "Optional, number of longitude points for the gaussian grid (recommend 2xLAT_S, default is for model to determine gaussian grid"
   "FHSTOCH", "Optional, forecast hour to write out random pattern in order to restart the pattern for a different forecast (used in DA), file is stoch_out.F<HHH>"
   "STOCHINI", "Optional, set to true if wanting to read in a previous random pattern (input file needs to be named stoch_ini)."

SPPT options 
""""""""""""

.. csv-table::
   :header: "Option", "Description"
   :widths: 30, 50

   "DO_SPPT", "logical to tell parent atmospheric model to use SPPT"
   "SPPT", "Amplitudes of random patterns (0.8,0.4,0.2,0.08,0.04) *"
   "SPPT_TAU", "Decorrelation timescales in seconds (21600,1.728E5,6.912E5,7.776E6,3.1536E7) *"
   "SPPT_LSCALE", "Decorrelation spatial scales in meters  (250.E3,1000.E3,2000.E3,2000.E3,2000.E3) *"
   "SPPT_LOGIT", "Should be true to limit the SPPT perturbations between 0 and 2.  Otherwise model will crash."
   "ISEED_SPPT", "Seeds for setting the random number sequence (ignored if stochini is true)"
   "SPPT_SIGTOP1", "lower sigma level to taper perturbations to zero (default is 0.1)"
   "SPPT_SIGTOP2", "upper sigma level to taper perturbations to zero (0.025)"
   "SPPT_SFCLIMIT", ".T.=tapers the SPPT perturbations to zero at model’s lowest level (helps reduce model crashes)"
   "SPPTINT", "Optional, interval in seconds to update random pattern.  Perturbations still get applied every time-step"
   "USE_ZMTNBLCK", ".T.=do not apply perturbations below the dividing streamline that is diagnosed by the gravity wave drag, mountain blocking scheme"

``*``  **SPPT** uses 5 different patterns of varying time/length scales that are added together before being passed to physics

SHUM options 
""""""""""""

.. csv-table::
   :header: "Option", "Description"
   :widths: 20, 50

   "DO_SHUM", "logical to tell parent atmospheric model to use SHUM"
   "SHUM", "Amplitudes of random patterns (0.004)"
   "SHUM_TAU", "Decorrelation timescales in seconds (21600)"
   "SHUM_LSCALE", "ecorrelation spatial scales in meters (250000)"
   "SHUM_SIGEFOLD", "e-folding lengthscale (in units of sigma) of specific humidity perturbations, default is 0.2)"
   "SHUMINT", "Optional, interval in seconds to update random pattern.  Perturbations still get applied every time-step"
   "ISEED_SHUM", "Seeds for setting the random number sequence (ignored if stochini is true)."

SKEB options
""""""""""""

.. csv-table::
   :header: "Option", "Description"
   :widths: 30, 50

   "DO_SKEB", "logical to tell parent atmospheric model to use SKEB"
   "SKEB", "Amplitudes of random patterns (0.5)"
   "SKEB_TAU", "Decorrelation timescales in seconds (21600)"
   "SKEB_LSCALE", "Decorrelation spatial scales in meters  (250)"
   "ISEED_SKEB", "Seeds for setting the random number sequence (ignored if stochini is true)."
   "SKEBNORM", "0-random pattern is stream function, 1-pattern is K.E. norm, 2-pattern is vorticity (default is 0)"
   "SKEB_VARSPECT_OPT", "0-gaussian (default), 1-power law (not tested)"
   "SKEB_NPASS", "number of passes of the del2 smoothing for the dissipation estimate (default is 11, minimum is 3)"
   "SKEB_VDOF", "the number of degrees of freedom in the vertical for the SKEB random pattern (default is 5)"
   "SKEB_SIGTOP1", "lower sigma level to taper perturbations to zero (default is 0.1)"
   "SKEB_SIGTOP2", "upper sigma level to taper perturbations to zero (0.025)"
   "SKEBINT", "Optional, interval in seconds to update random pattern.  Perturbations still get applied every time-step"

