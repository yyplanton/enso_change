# Correcting underestimated internal variability fails to reconcile models with observed Pacific SST gradient strengthening

Submitted to [Geophysical Research Letters](https://agupubs.onlinelibrary.wiley.com/journal/19448007?gad_campaignid=20631722826&gad_source=1&gclid=EAIaIQobChMIwtTU0eyXkAMVAjStBh0bXDBkEAAYASAAEgIDhPD_BwE&utm_campaign=R3MR425&utm_content=PMAX&utm_medium=cpc&utm_source=google)

## Getting started
#### Install
This method uses [pip](https://pip.pypa.io/en/stable/installation/)  
Download the repository: `git clone <repo>`  
Enter the directory: `cd enso_and_trends`  
If you haven't already installed virtualenv: `pip install virtualenv`  
Create your new environment (called 'enso_and_trends'): `virtualenv enso_and_trends`  
Activate your new environment: `source enso_and_trends/bin/activate`  
Install the requirements in the current environment: `pip install -r requirements.txt`
#### Data
Data for all figures but Figure 1 is already available with the repository.  
To access data for Figure 1 (i.e., maps of SST trends from 64 CMIP6 models), download from: [Harvard Dataverse](https://doi.org/10.7910/DVN/FOVDTN)   
## Paper
#### Links
URL or DOI: [preprint](https://doi.org/10.22541/au.176125924.49048026/v1)
#### Abstract
The equatorial Pacific zonal sea surface temperature (SST) gradient has strengthened since 1980, yet fewer than 1% of CMIP6 simulations reproduce this trend. We test whether underestimated internal variability explains this mismatch. Extreme El Niño events enhance interdecadal variability of the gradient, but CMIP6 models simulate them too infrequently, producing low Niño3 SST skewness that correlates (~0.65) with interdecadal gradient variability. Applying this observational constraint to statistically correct modeled interdecadal gradient variability increases the fraction of simulations exceeding the observed trend to only ~6%, indicating that enhancing variability alone cannot resolve the discrepancy. One-third of models with large ensembles display a weak but forced recent strengthening that later reverses. Beyond +1.7°C of global warming, over 95% of simulations project a ~15% weakening. Together, these results suggest that the recent strengthening includes a substantial externally forced transient component poorly simulated by models, while the projected long-term weakening is robust.
#### Figure 1
Observed and CMIP6 SST trends over 1980–2014, with the tropical mean (20°S–20°N) removed to highlight spatial patterns. a) Observations, b) CMIP6 multimodel mean (MMM), mean of members in the c) bottom and d) top 10% of the zonal SST gradient trend distribution, e) ensemble mean of 12 ensembles with a forced strengthening component, f) Meridional averages (5°S–5°N) of panels a–e.
#### Figure 2
Influence of extreme El Niño events on interdecadal zonal SST gradient trends. a) Niño3 SST anomalies and b) zonal SST gradient, with red shading marking extreme El Niño events (April to March 1982/83, 1997/98, 2015/16; Santoso et al., 2017). In (a), g1 indicates Niño3 skewness over 1980–2014 using all data (black) and excluding extremes (red). In (b), dashed black (all years) and solid red (excluding extremes) lines show 1980–2014 gradient trends, with slopes given below. Δr(σ[s]) is the percentage reduction in the standard deviation in all possible 35-year gradient trends over 1980–2024, when extremes are excluded.
#### Figure 3
Observational constraint linking extreme El Niño events to the amplitude of interdecadal zonal SST gradient trends. a) Ensemble standard deviation of 35-year gradient trends versus ensemble mean Niño3 SST anomaly skewness in CMIP6 (36 historical ensembles for 1980–2014 in green, 12 historical ensembles with a forced strengthening component during 1980–2014 in dark blue, and 64 piControl ensembles in purple). Solid lines show linear regressions; correlations and 2.5-97.5% confidence intervals are given at bottom right. The vertical line marks the observed 1980–2014 Niño3 skewness. The dashed horizontal line shows the observed standard deviation of 35-year gradient trends (computed from 5-year-staggered 35-year windows over 1870–2024) and the grey shading marks the observational range. b) Corresponding historical (green) and historical with a forced strengthening component (dark blue) CMIP6 trend distributions (thick line: median; dashed: 25th–75th percentiles; whiskers: 1st–99th). The horizontal line marks the observed 1980–2014 gradient strengthening. “Corrected” distributions are broadened using the regression in (a) and the model–observations skewness mismatch (supplementary material Text S2). Percentages indicate the fraction of each distribution exceeding the observed value.
#### Figure 4
Global warming temperature level at which zonal SST gradient weakening emerges in CMIP6. a) Percentage of members (averaged across individual model ensembles) with a gradient below the piControl mean as a function of GMST anomaly. b) Change in the zonal gradient (%) relative to the piControl mean. Lines (shading) show the median (25th–75th percentiles) for all CMIP6 models (green), the bottom 10% (dashed red), top 10% (dot-dashed blue) ensemble members ranked by their 1980–2014 historical trends, and 12 historical ensembles with a forced strengthening component during 1980–2014 in (dark blue). c) Observed (black, relative to 1850-1900) and modeled (green for historical, other colors for various SSPs, relative to piControl) global warming level. Lines (shading) show the mean (25th–75th percentiles; 95% confidence interval for observations). Vertical lines on a,b mark the ~0.6 °C global warming level at the end of the historical period (1980–2014), and the +1.7°C warming level.
#### Definitions
- Extreme El Niño events:
  * observed after 1980 are defined as 1982/83, 1997/98, and 2015/16 (Santoso et al.,
[2017](https://doi.org/10.1002/2017RG000560)).
- Interannual anomalies:
  * monthly mean seasonal cycle (computed over a given epoch) removed
- Regions:
  * eastern equatorial Pacific (E): [5S-5N ; 80-180W]
  * Niño3 (N3): [5S-5N ; 90-150W]
  * tropic: [20S-20N ; 0-360E]
  * western equatorial Pacific (W): [5S-5N ; 110-180E]
  * zonal Pacific gradient (&Delta;<sub>W-E</sub>): values averaged over western equatorial Pacific minus values
  averaged over eastern equatorial Pacific
- Statistics:
  * Skewness (Cramér, [1946](https://doi.org/10.1515/9781400883868)): 
$g_1 = {1/n \sum_{i=1}^n (x_i - \bar{x}) \over [1/n \sum_{i=1}^n (x_i - \bar{x})^2]^{2/3}} $
  * Trend: slope (&alpha;) of the linear regression computed along the time dimension $y = \alpha.x +\beta$
- Variables:
  * global surface air temperature (GSAT) from [HadCRUT5](https://www.metoffice.gov.uk/hadobs/hadcrut5/)
  * sea surface temperature (SST) from [COBE-SST 2](https://psl.noaa.gov/data/gridded/data.cobe2.html), 
  [ERSSTv5](https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html) and 
  [HadISST](https://www.metoffice.gov.uk/hadobs/hadisst/)
  * surface temperature (TS) from 64 [CMIP6](https://esgf‐node.llnl.gov/projects/esgf‐llnl/) model ensembles
