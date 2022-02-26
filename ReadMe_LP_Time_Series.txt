Geochemical X ray fluorescence log ratio time series data for two sediment cores, LP08 and LP16, extracted from Lago Pato, Torres del Paine, Southern Chile 

Abstract
The dataset comprises of  X ray fluorescence log ratio time series data for two sediment cores from Lago Pato, a small lake basin at S51°18.020’, W72°40.716’ and ~ 33 m a.s.l., which is topographically separated from Lago del Toro in Torres del Paine (TdP). The data are used to constrain glacier dynamics and lake level change in the TdP and Última Esperanza region over the last ~30,000 cal a BP (30 ka). LP08 was extracted from the current depocentre in November 2007 to March 2008. LP16 was extracted the terrestrial shoreline in November 2015.

Funding source
This project was funded by the Natural Environment Research Council (NERC) through the British Antarctic Survey (BAS) and an UGent BOF bilateral collaboration project. RMcC was supported by Programa Regional R17A10002 and R20F0002 (PATSER) ANID. We gratefully acknowledge the University of Magallanes (UMAG) and the University of Santiago (Carolina Diaz) for assistance with fieldwork; the NERC/SUERC AMS Radiocarbon Facility for providing initial range-finder radiocarbon dates; the NERC Isotope Geosciences Laboratory (NIGL, now National Environmental Isotope Facility, NEIF, at the British Geological Survey) and Melanie Lang for stable carbon isotope analysis; Aberystwyth University (David Kelly), Durham University (Neil Tunstall and Christopher Longley) and Edinburgh University (Chris Hayward) for use of their core scanning and microprobe facilities and technical support.

Keywords
Last Glacial Maximum (LGM), palaeoclimate, palaeolimnology, glaciation, lake level changes, Patagonia, Southern Hemisphere Westerly Winds.

Personnel
Data collectors (ORCID code)
Stephen J. Roberts1 (0000-0001-5542-3703) – sediment core extraction, chronology, geochemistry
Robert D. McCulloch2 (0000-0001-5542-3703) – chronology
Joseph F. Emmings3 – geochemistry
Sarah J. Davies4 – geochemistry

Data analysts
Stephen J. Roberts1 – all data
Robert D. McCulloch2 – chronology
Joseph F. Emmings3 – geochemistry
Sarah J. Davies4 – geochemistry

Affiliations
1British Antarctic Survey, Natural Environment Research Council, High Cross, Madingley Road, Cambridge, CB3 0ET, UK.
2Centro de Investigación en Ecosistemas de la Patagonia (CIEP), Coyhaique, Chile.
3British Geological Survey, Keyworth, Nottingham, NG12 5GG, UK.
4Department of Geography and Earth Sciences, Aberystwyth University, Aberystwyth, SY23 3DB, UK.

Methodology

Chronology

A chronology for each record was established using Accelerator Mass Spectrometry (AMS) radiocarbon dating of 21 samples from the LP08 record and 15 samples from the LP16 record. Calibration of radiocarbon ages was undertaken in OXCAL v.4.4 using the SHCal20.14C Southern Hemisphere atmosphere calibration curve (SH20). Radiocarbon ages are reported as conventional radiocarbon years BP (14C years BP) ±1σ and calibrated ages as 2σ (95.4%) ranges, median and mean calendar years BP (cal a BP and cal ka BP, relative to 1950 CE), rounded to the nearest ten years. Age-depth models were developed using Bayesian age-depth modelling software (rBACON v.2.5). Modelled age data mean ages produced by the SH20M1H (Southern Hemisphere, SHcal20, radiocarbon calibration curve) in rBACON, where M1 indicates Model 1 and H indicates the inclusion of a hiatus in the model.

Geochemistry
Sediment cores were collected using a UWITEC-gravity corer, Livingston piston corer and a Russian corer from the deepest point (~3.5 m of water depth) in Lago Pato:
-	LP08 record from  S51° 18’01.2’’, W72° 40’43.0’’, 32 m a.s.l. is 600 cm long
-	LP16 record from S51°18’11.3’’, W72° 40’53.7’’, 33–34 m a.s.l. is 295 cm long

Contiguous downcore wet-sediment Energy Dispersive Spectrometry (EDS) X-ray fluorescence core scanning (XRF-CS) data was collected using an ITRAX XRF core scanner at Aberystwyth University fitted with a Molybdenum (Mo) anode X-ray tube (settings: 30 kV, 50 mA, count time 10 seconds, at 2 mm contiguous intervals and for LP08 Unit 6 (equivalent to mean ± 2-sigma: 4.5±7.0 years), at 200 μm intervals for LP08 Unit 1 (1.3±4.2 years),  with LP08 basal Unit 1 scanned at 100 μm, and at 500 μm for LP16 Units 2-6 (9.6±17.4 years) and at 200 μm for LP16 Unit 1 200 μm (1.1±1.6 years). Data from finely laminated glaciolacustrine sediments in Units 1-2 were measured at or smoothed to 200 μm (from 100 μm interval data) before analysis.


Time series analysis 
Log-n element/Ti ratio XRF-CS Z-scores were used for time series analysis (Fast Fourier Transform, FFT, periodograms, Lomb-Scargle Power Spectrum, Wavelet Power Spectrum, Peak Identification) in MATLAB. Equally spaced (10-year and 100-year) time-intervals were generated using a Piecewise Cubic Hermite Interpolated Polynomial (PCHIP) function, which avoids spline artefacts and preserve the shape of the original XRF-CS data series (Grinsted et al., 2004; Trauth, 2015). Time series data were detrended (polynomial linear best fit) to remove the long-term linear trend. Second order polynomial Locally Weighted Scatterplot Smoothing (LOESS) 100-year smoothing (0.1 sampling interval with outliers removed) was also used to compare datasets to published data.

Instrumentation
Data were analysed in MATLAB v. R2021a, R v. 4.1.0/Rstudio v. 1.4.171,  using the R packages Vegan, Rioja, Tidyverse, ggplot2, Ggally v. 2.1.2. Code is available from: https://github.com/stever60/Lago_Pato 

Quality
NA.

Related URLs
Code, input data and output summary figures are available on: www.github.com/stever60 

Temporal coverage
The LP08 and LP16 records cover the time period between 0–30,000 years

Spatial coverage
Sediment core data are from two sites:
-	LP08 record from  S51° 18’01.2’’, W72° 40’43.0’’, 32 m a.s.l. is 600 cm long
-	LP16 record from S51°18’11.3’’, W72° 40’53.7’’, 33–34 m a.s.l. is 295 cm long

Resolution
XRF-CS data from finely laminated glaciolacustrine sediments in Units 1-2 were smoothed to 200 μm, other data to 2mm, and equal spaced time-intervals (10-years and 100-years) for use in time series analysis were generated. Time series datasets were centred, standardised (Z-scores), detrended and interpolated using a 10-year or 100-year Piecewise Cubic Hermite Interpolating Polynomial (PCHIP) function, which avoids spline artefacts and preserve the shape of the original XRF-CS data series.

Location
Lago Pato, Torres del Paine National Park, Patagonia, Chile

References
Roberts SJ, McCulloch RDM, Emmings J, Davies SJ, Van Nieuwenhuyze W, Sterken M, et al. (2022) Late glacial and Holocene palaeolake history of the Última Esperanza region of Southern Patagonia. Frontiers in Earth Science.

Grinsted, A., Moore, J.C., and Jevrejeva, S. (2004). Application of the cross wavelet transform and wavelet coherence to geophysical time series. Nonlinear Processes in Geophysics 11, 561-566.

Trauth, M.H. (2015). "Time-Series Analysis," in MATLAB® Recipes for Earth Sciences, ed. M. H. Trauth.  (Berlin, Heidelberg: Springer Berlin Heidelberg), 151-213.

Data structure and data format
The LP08 and LP16 datasets are arranged as follows:

The input output csv files are arranged into subfolders as follows. These subfolders represent different time periods in the LP08 and LP16 records. Changing the time period investigated alters the Z-scores produced. 

LP08_5ka - covering the last 5000 years
LP08_8ka - covering the last 8000 years
LP08_10ka - covering the last 10,000 years
LP08_Unit1 - covering 21,180-29,780 cal a BP
LP08_Unit1_basal - covering 26,490-29,780 cal a BP

LP16_11ka - covering the last 11,000 years
LP16_14ka - covering the last 14,000 years
LP16_Unit1 - covering 20,400-27,550 cal a BP

LP08_LP16_10ka - comparison of dataset pairs for LP08 and LP16 records covering the last 10,000 years
LP08_LP16_21_27ka - comparison of dataset pairs for LP08 and LP16 records covering the 21-27 ka cal BP
LP08_LP16_Fig9_LnFe_Mn – Fe/Mn natural log ratio time series data used in the Frontiers paper for LP008 and LP16 records

The following log ratio dataset pairs were run:
1_Fe_Mn_&_Mn_Ti
2_Fe_Mn_&_Br_Ti
3_Mn_Ti_&_Br_Ti
4_Fe_Mn_&_Inc_Coh
5_Mn_Ti_&_Ca_Ti
6_Br_Ti_&_Inc_Coh

Input datafiles:

The filename within each folder matches the folder name and indicates the time period covered by the time series. For example:

LP08_Unit1_basal_inputs.csv

Output datafiles (e.g., LP08_Unit1_basal_ouput.csv)

The filename within each folder matches the folder name and indicates the time period covered by the time series. 

The key for the columns in each dataset pair (e.g., 1_Fe_Mn_&_Mn_Ti)  are as follows (in order from left to right. 

% As measured data	
s1x, s1xn	series 1x time scale (cal a BP)
s2x, s2xn	series 2x time scale (cal a BP)
series1_input	series 1y input data (log ratios)
series2_input	series 2y input data (log ratios)
s1n	 	series 1y normalised (Z-scores)
s2n	 	series 2y normalised (Z-scores)
series1ydl	series 1y detrended - linear best fit
series1ydlnorm	series 1y standardised (Z-scores) detrended (linear best fit)
series2ydl	series 2y detrended - linear best fit
series2ydlnorm	series 2y standardised (Z-scores) detrended (linear best fit)
	
% Interpolated data	
sx1inter	series  	1x time scale interpolated (10 or 100 year intervals)
sy1inter	series  	1y interpolated (10 or 100 year intervals)
sy1norminter	series 1y standardised (Z-scores) & interpolated (10 or 100 year intervals)
sx2inter	series  	2x time scale interpolated (10 or 100 year intervals)
sy2inter	series  	2y interpolated (10 or 100 year intervals)
sy2norminter	series 2y standardised (Z-scores) & interpolated (10 or 100 year intervals)
	
% Detrended & interpolated data (s = reshaped array into wide format)
x2 or x2s	series 1x time scale interpolated (10 or 100 year intervals)
x4 or x4s	series 2x time scale interpolated (10 or 100 year intervals)
y2 or y2s 	series 1y detrended and factor interpolated data - subtracted from mean
y2l or y2sl 	series 1y detrended and factor interpolated data - subtracted from linear
y4 or y4s 	detrended and factor interpolated data series 2 - subtracted from mean
y4l or y4sl 	detrended and factor interpolated data series 2 - subtracted from linear
	
% Standardized, Detrended & interpolated data	
y2ln or y2sln 	series 1y standardised and detrended data - subtracted from linear best fit
y4ln or y4sln 	series 2y standardised (Z-scores) and detrended and factor interpolated data - subtracted from linear
	
% Standardized, Detrended & interpolated periodicity data	
s1period	 periodicity of series 1y standardized (Z-scores) and detrended (subtracted from linear best fit) interpolated data (i.e., y2ln or y2sln)
s2period	 periodicity of series 2y standardized (Z-scores) and detrended (subtracted from linear best fit) interpolated data (i.e., y4ln or y4sln)
msc	 	Magnitude squared coherence (MSC) of s1period and s2period for the time series
frcspect	 	Frequency cross spectrum (FCS) of s1period and s2period for the time series

Access constraints
None

Use constraints
NERC-funded data, the Open Government Licence applies.


