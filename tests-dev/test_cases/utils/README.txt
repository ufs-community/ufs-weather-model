2020 CAPE bash plotting script developed by Ratko Vasic.

Baroclinic wave bash plotting script developed by Ratko Vasic,
adapted from an original grads script & configuration files provided by Xiaqiong Zhou.



To generate plots, copy the relevant .sh script into your UFS-WM test run directory 
(i.e., `run_dir/<case_name>_<compiler>` from within ufs-weather-model/tests-dev). 
That is, to generate plots for a 2020 CAPE test, copy `plot_cape.sh` into your
2020_CAPE_<compiler> run directory; to generate plots for the baroclinic wave case,
copy `plot_bcw.sh` into your <baroclinic_wave_<compiler> run directory.

To run the scripts, invoke from the command line via `./plot_<bcw/cape>.sh`. 

Users can adjust the experiment name, standard pressure level, and forecast hour in the baroclinic wave script. 
For the CAPE script, users can choose between global/regional(CONUS) domain, as well as which 
analysis file resolution to generate plots for.
