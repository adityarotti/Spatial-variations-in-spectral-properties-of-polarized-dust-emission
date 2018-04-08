# Spatial-variations-in-spectral-properties-of-polarized-dust-emission

To do list:
  - Simulate a map with a fixed spectral slope and some amplitude. 
  - Modulate this map to generate amplitude variations across the maps.(Maybe this is not necessary) 
  - Add Planck noise and beam smoothing. 
  - Now use this synthetic data to recover a map of amplitude and slope as a base model following exactly the same analysis procedure as that on data. 
  - Compare the results (specifically the statistics of the inferred amplitude and slope) of this base model to ones derived from Planck.
  
This exercise will give us a sense of how much fluctuation is expected from Planck noise model alone. This exercise may even be repeated for a few 10 simulations to generate a reasonable statistic. 

In a sense we are claiming that there are some regions of the sky which do not follow this base model. 


Questions:

 -  How do you associate ell to binned spectra? Right now we are using the mean of the multipoles that go into the bin. But clearly the fitted slope depends on this. Is there a well defined way to do it? Or do we need to take into account the width of the bin into the fit ? 



Have run some preliminary analysis to visualize the parameter covariance. More such plots can be seen @ https://github.com/adityarotti/Spatial-variations-in-spectral-properties-of-polarized-dust-emission/blob/master/Analysis/analysis_april2018/results/disc_analysis/gal80/obs/nside4_lmin60_lmax600_dl40_mlmax256_discsize20.0_apow2/parameter_plots/

<img src="https://github.com/adityarotti/Spatial-variations-in-spectral-properties-of-polarized-dust-emission/blob/master/Analysis/analysis_april2018/results/disc_analysis/gal80/obs/nside4_lmin60_lmax600_dl40_mlmax256_discsize20.0_apow2/parameter_plots/amp_slope_ipix0.jpg" alt="E-mode map 353 GHz" width="48%"></img>
