Matlab Time Series Analysis 

Matlab code, input/output data and summary figures can be found here: https://github.com/stever60/Lago_Pato 

Age data is in cal years before present (cal a BP), where present is 1950 CE. Ages data are modelled and based on the mean SH20M1H (Southern Hemisphere, SHcal20, radiocarbon calibration curve) produced by the M1 (Model1) Bayesian depth model in rBACON. H indicates the inclusion of a 10,000 year hiatus at 470 cm depth in the model (SH20M1H). 

Input and output datafiles data are arranged into folders as follows. These represent different time periods in the LP08 and LP16 records. Changing the time period investigated alters the Z-scores produced. 

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

Frontiers - time series data presented in the Frontiers paper for LP008 and LP16 records

The following log ratio dataset pairs were run and arranged into folders as follows:
1_Fe_Mn_&_Mn_Ti
2_Fe_Mn_&_Br_Ti
3_Mn_Ti_&_Br_Ti
4_Fe_Mn_&_Inc_Coh
5_Mn_Ti_&_Ca_Ti
6_Br_Ti_&_Inc_Coh
Frontiers

Input datafiles:

The filename within each folder matches the folder name and indicates the time period covered by the time series. For example:

LP08_Unit1_basal_inputs.csv


Output datafiles (e.g., LP08_Unit1_basal_ouput.csv)

The filename within each folder matches the folder name and indicates the time period covered by the time series. The key for the columns in each dataset pair (e.g., 1_Fe_Mn_&_Mn_Ti)  are as follows (in order from left to right. 

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
s1period	 	periodicity of series 1y standardized (Z-scores) and detrended (subtracted from linear best fit) interpolated data (i.e., y2ln or y2sln)
s2period	 	periodicity of series 2y standardized (Z-scores) and detrended (subtracted from linear best fit) interpolated data (i.e., y4ln or y4sln)
msc	 	Magnitude squared coherence (MSC) of s1period and s2period for the time series
frcspect	 	Frequency cross spectrum (FCS) of s1period and s2period for the time series












