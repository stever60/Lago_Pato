%% Clear everything and start again
clc
clear all
close all
 
%%  Add paths for mac
addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/Wavelet/grinsted-wavelet-coherence-d987ea4');
addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/Wavelet/grinsted-wavelet-coherence-d987ea4/data');
addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis');
addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Data')
addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Data/LP16')
addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Code')
addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Code/LP16')
addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/fLOESS')

%% Useful code
%remove NaNs code
%out = A(all(~isnan(A),2),:); % for nan - rows
%out = A(:,all(~isnan(A)));   % for nan - columns

%% PART 1: IMPORT DATA, CHECK, CREATE Z-scores 

% Load datafiles into matlab
% Run one paired dataset at a time & uncomment/comment below out
% Imported data have imcomplete 10-yr top and bottom lines removed prior to
% import to prevent extrapolation beyond dataset min and max limits

%% DATASET 1

% Series 1 & 2 and remove NaN rows
%load LP16_U1_LnFe_Mn.txt
%series1 = LP16_U1_LnFe_Mn(all(~isnan(LP16_U1_LnFe_Mn),2),:);
%load LP16_U1_LnMn_Ti.txt
%series2 = LP16_U1_LnMn_Ti(all(~isnan(LP16_U1_LnMn_Ti),2),:);

% Create series labels and x, y axis labels for all plots
%seriesname1={'LP16 U1 Ln(Fe/Mn)'};
%ylbl1 = {'LP16 U1: Ln(Fe/Mn)'};
%xlbl1 = {'Age (cal yr BP)'};
%seriesname2={'LP16 U1 Ln(Mn/Ti)'};
%ylbl2 = {'LP16 U1: Ln(Mn/Ti)'};
%xlbl2 = {'Age (cal yr BP)'};
%ylblspecial1 = {'LP16 U1: Ln(Fe/Mn) (blue), LP16 U1: Ln(Mn/Ti) (orange)'};

%% DATASET 2
% Series 1 and remove NaN row
%load LP16_U1_LnFe_Mn.txt
%series1 = LP16_U1_LnFe_Mn(all(~isnan(LP16_U1_LnFe_Mn),2),:);
%load LP16_U1_LnBr_Ti.txt
%series2 = LP16_U1_LnBr_Ti(all(~isnan(LP16_U1_LnBr_Ti),2),:);

% Create series labels and x, y axis labels for all plots
%seriesname1={'LP16 U1 Ln(Fe/Mn)'};
%ylbl1 = {'LP16 U1: Ln(Fe/Mn)'};
%xlbl1 = {'Age (cal yr BP)'};
%seriesname2={'LP16 U1 Ln(Br/Ti)'};
%ylbl2 = {'LP16 U1: Ln(Br/Ti)'};
%xlbl2 = {'Age (cal yr BP)'};
%ylblspecial1 = {'LP16 U1a: Ln(Fe/Mn) (blue), LP16 U1: Ln(Br/Ti) (orange)'};

%% DATASET 3
% Series 1
%load LP16_U1_LnMn_Ti.txt
%series1 = LP16_U1_LnMn_Ti(all(~isnan(LP16_U1_LnMn_Ti),2),:);
%load LP16_U1_LnBr_Ti.txt
%series2 = LP16_U1_LnBr_Ti(all(~isnan(LP16_U1_LnBr_Ti),2),:);

% Create series labels and x, y axis labels for all plots
%seriesname1={'LP16 U1 Ln(Mn/Ti)'};
%ylbl1 = {'LP16 U1: Ln(Mn/Ti)'};
%xlbl1 = {'Age (cal yr BP)'};
%seriesname2={'LP16 U1 Ln(Br/Ti)'};
%ylbl2 = {'LP16 U1: Ln(Br/Ti)'};
%xlbl2 = {'Age (cal yr BP)'};
%ylblspecial1 = {'LP16 U1: Ln(Mn/Ti) (blue), LP16 U1: Ln(Br/Ti) (orange)'};

%% DATASET 4
% Series 1
%load LP16_U1_LnFe_Mn.txt
%series1 = LP16_U1_LnFe_Mn(all(~isnan(LP16_U1_LnFe_Mn),2),:);
%load LP16_U1_LnInc_Coh.txt
%series2 = LP16_U1_LnInc_Coh(all(~isnan(LP16_U1_LnInc_Coh),2),:);

% Create series labels and x, y axis labels for all plots
%seriesname1={'LP16 U1 Ln(Fe/Mn)'};
%ylbl1 = {'LP16 U1: Ln(Fe/Mn)'};
%xlbl1 = {'Age (cal yr BP)'};
%seriesname2={'LP16 U1 Ln(Inc/Coh)'};
%ylbl2 = {'LP16 U1: Ln(Inc/Coh)'};
%xlbl2 = {'Age (cal yr BP)'};
%ylblspecial1 = {'LP16 U1: Ln(Fe/Mn) (blue), LP16 U1: Ln(Inc/Coh) (orange)'};

%% DATASET 5
% Series 1
%load LP16_U1_LnMn_Ti.txt
%series1 = LP16_U1_LnMn_Ti(all(~isnan(LP16_U1_LnMn_Ti),2),:);
%load LP16_U1_LnCa_Ti.txt
%series2 = LP16_U1_LnCa_Ti(all(~isnan(LP16_U1_LnCa_Ti),2),:);

% Create series labels and x, y axis labels for all plots
%seriesname1={'LP16 U1 Ln(Mn/Ti)'};
%ylbl1 = {'LP16 U1: Ln(Mn/Ti)'};
%xlbl1 = {'Age (cal yr BP)'};
%seriesname2={'LP16 U1 Ln(Ca/Ti)'};
%ylbl2 = {'LP16 U1: Ln(Ca/Ti)'};
%xlbl2 = {'Age (cal yr BP)'};
%ylblspecial1 = {'LP16 U1: Ln(Mn/Ti) (blue), LP16 U1: Ln(Ca/Ti) (orange)'};

%% DATASET 6
% Series 1
load LP16_U1_LnBr_Ti.txt
series1 = LP16_U1_LnBr_Ti(all(~isnan(LP16_U1_LnBr_Ti),2),:);
load LP16_U1_LnInc_Coh.txt
series2 = LP16_U1_LnInc_Coh(all(~isnan(LP16_U1_LnInc_Coh),2),:);

% Create series labels and x, y axis labels for all plots
seriesname1={'LP16 U1 Ln(Br/Ti)'};
ylbl1 = {'LP16 U1: Ln(Br/Ti)'};
xlbl1 = {'Age (cal yr BP)'};
seriesname2={'LP16 U1 Ln(Inc/Coh)'};
ylbl2 = {'LP16 U1: Ln(Inc/Coh)'};
xlbl2 = {'Age (cal yr BP)'};
ylblspecial1 = {'LP16 U1: Ln(Fe/Mn) (blue), LP16 U1: Ln(Inc/Coh) (orange)'};

%% Set up two series to run 

% Create common file name and take natural log of the data if needed
% y = log(x) returns the natural logarithm ln(x) of each element in the array

%Series 1
format long;
series1x = (series1(:,1));
series1y = (series1(:,2));
%series1y = log(Series1(:,2));
x1 = series1x;
y01 = series1y; %original data series 1

%check output - if needed - ; stops the output displaying in command window
series1x
series1y

% Series 2
format long;
series2x = (series2(:,1));
series2y = (series2(:,2));
%series2ly = log(series2(:,2));
x3 = series2x;
y02 = series2y; %original data series 2
series2x
series2y

% check reshape value here for interpolation
rs = 715

%% Standardize and centre y (Z-scores) 
% subtracting from mean and dividing by stdev for series 1 and 2
series1ynorm = (series1y-mean(series1y))/std(series1y)
series2ynorm = (series2y-mean(series2y))/std(series2y)

%write to file 
s1norm = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/s1n.csv';
xlswrite(s1norm,series1ynorm);
s1xnorm = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/s1xn.csv';
xlswrite(s1xnorm,series1x);
s2norm = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/s2n.csv';
xlswrite(s2norm,series2ynorm);
s2xnorm = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/s2xn.csv';
xlswrite(s2xnorm,series2x);

%% Plot to check input series look OK
figure(1)
ax(1) = subplot(2,1,1);
plot(series1x,series1y,series1x,series1ynorm);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
title(['\fontsize{12}','Time series plot for ' seriesname1{1}]);
xlim([min(series1x) max(series1x)]), hold on
hold off
xlim auto;
ylim auto;
legend('Measured', 'Z-scores');
legend('Location','eastoutside')
legend ('boxon');

ax(2) = subplot(2,1,2);
plot(series2x,series2y,series2x,series2ynorm);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
title(['\fontsize{12}','Time series plot for ' seriesname2{1}]);
xlim([min(series2x) max(series2x)]);
xlim auto;
ylim auto;
legend('Measured', 'Z-scores');
legend('Location','eastoutside')
legend ('boxon');

%link parts and zoom in to show the changes
linkaxes(ax(1:2),'x');
%xlim([min(series1x) max(series1x)])
xlim ([min(-100) max(20000)])
xlim auto
ylim auto

h1 = figure(1)
set(h1, 'Units', 'centimeters','PaperSize', [20,20]);
print(h1, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig1', '-dpdf')

%% PART 2: INTERPOLATION OF TIME SERIES 

% Establish min and max limits for x-axis (Years) and interval spacing 

%Series 1
format long;
yrmin1 = round(min(series1x)) %rounded to the nearest year
yrmax1 = round(max(series1x)) %rounded to the nearest year
intv1 = diff(series1x)
intv1sm = smoothdata(intv1)
meanyr1 = mean(intv1)
meanyr1sm = smoothdata(intv1sm)
stdevyr1 = 2*std(intv1)

%Series 2
format long;
yrmin2 = round(min(series2x)) 
yrmax2 = round(max(series2x))
intv2 = diff(series2x)
intv2sm = smoothdata(intv2)
meanyr2 = mean(intv2)
meanyr2sm = smoothdata(intv2sm)
stdevyr2 = 2*std(intv2)

%Look at mean interval for series 1 and series 2 in years
yrmin1
yrmax1
yrmin2
yrmax2
meanyr1
meanyr2
stdevyr1
stdevyr2

% interpolation produces an equally-spaced time series equal to the mean time 
% series interval of the dataset
% For LP08_Unit 1 data, this case can be as low as 1 
% because mean interval is 0.6643 years per sample point
% but have set the mean interval to 10 years to remove low frequency noise

%% Calculate smoothing factors

% convert the time interval into an integer by rounding to nearest year  
meanyr1sm = mean(intv1sm)
meanyr2sm = mean(intv2sm)
meanyr1round = round(meanyr1)
meanyr2round = round(meanyr2)

% Smooth to nearest 10 years
% 10-yr dataset - need to enter factor = 10 as data input here for 10 yr interpolation
factor1 = round(meanyr1*(10/meanyr1)) %100
factor2 = round(meanyr2*(10/meanyr2)) %100
meanyr1x = factor1
meanyr2x = factor2

%1-yr dataset to make two vectors of the same size
factor3 = round(meanyr1*(1/meanyr1))
factor4 = round(meanyr2*(1/meanyr2))
meanyr3x = factor3
meanyr4x = factor4

%100-yr dataset to make two vectors of the same size
factor5 = round(meanyr1*(100/meanyr1))
factor6 = round(meanyr2*(100/meanyr2))
meanyr5x = factor5
meanyr6x = factor6

%% Plot time interval summary and add values

figure(2)
ax(1) = subplot(2,1,1);
plot(intv1)
hold on
plot(intv1sm)
title(['\fontsize{12}','Time interval plot for ' seriesname1{1}]);
xlabel('Sample no.');
ylabel(ylbl1{1});
ylim auto %([min(0) max(20)]);
legend(['Mean interval (years±2s) = ',num2str(meanyr1),  '±',num2str(stdevyr1)],['Mov. Ave. mean = ',num2str(meanyr1sm),  ' years']);
legend('Location','southoutside')
legend ('boxon');
hold off

ax(2) = subplot(2,1,2);
plot(intv2)
hold on
plot(intv2sm)
title(['\fontsize{12}','Time interval plot for ' seriesname2{1}]);
xlabel('Sample no.');
ylabel(ylbl2{1});
ylim auto %([min(0) max(20)]);
legend(['Mean interval (years±2s) = ',num2str(meanyr2),  '±',num2str(stdevyr2)],['Mov. Ave. mean = ',num2str(meanyr2sm),  ' years']);
legend('Location','southoutside')
legend ('boxon');
hold off

%link parts and zoom in to show the changes
linkaxes(ax(1:2),'x');
%xlim([min(series1x) max(series1x)])
%xlim ([min(-100) max(20000)])
xlim auto
ylim auto

h2 = figure(2)
set(h2, 'Units', 'centimeters','PaperSize', [20,20]);
print(h2, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig2', '-dpdf')

%% Work out new x-series for interpolated data  

%define t1 and t2 as yrmin: mean yr interval: yr max series needed for interpolation
%to create evenly-spaced time series without increasing the number of
%datapoints
t1x = yrmin1:meanyr1x:yrmax1;
t2x = yrmin2:meanyr2x:yrmax2;

t1 = round(t1x,-1)
t2 = round(t2x,-1)

%1-yr PCHIP interpolated dataset 
t3x = yrmin1:meanyr3x:yrmax1;
t4x = yrmin2:meanyr4x:yrmax2;

t3 = round(t3x,-1)
t4 = round(t4x,-1)

% Value for T1 is used later on - don't delete
%t1 = yrmin1 : T1 : yrmax1;
%t2 = yrmin2 : T2 : yrmax2;
T1 = meanyr1;
T2 = meanyr2;
Fs1 = 1/T1; 
Fs2 = 1/T2;

min(t1)
max(t1)
t1;
t2;

%% Interpolate Series 1 and 2 data and standardized (not detrended version)
%PCHIP - mean
s1inter1 = interp1(series1x,series1y,t1,'pchip'); %t1 for adjusted time series = x2
s2inter1 = interp1(series2x,series2y,t2,'pchip'); %t1 for adjusted time series = x2

s1norminter = interp1(series1x,series1ynorm,t1,'pchip'); %t1 for adjusted time series = x2
s2norminter = interp1(series2x,series2ynorm,t2,'pchip'); %t1 for adjusted time series = x2

%reshape interpolated y series 1 and 2 into a vector with n rows x 1
%column from n columns x 1 row - to use later on for regression etc.

rs = 715

x1reshape = reshape(t1,[rs,1]);
x2reshape = reshape(t2,[rs,1]);

y1reshape = reshape(s1inter1,[rs,1]);
y2reshape = reshape(s2inter1,[rs,1]);
y1normreshape = reshape(s1norminter,[rs,1]);
y2normreshape = reshape(s2norminter,[rs,1]);

%write to file 
sx1inter = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/sx1inter.csv';
xlswrite(sx1inter,x1reshape);
sx2inter = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/sx2inter.csv';
xlswrite(sx2inter,x2reshape);

sy1inter = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/sy1inter.csv';
xlswrite(sy1inter,y1reshape);
sy2inter = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/sy2inter.csv';
xlswrite(sy2inter,y2reshape);

sy1norminter = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/sy1norminter.csv';
xlswrite(sy1norminter,y1normreshape);
sy2inter = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/sy2norminter.csv';
xlswrite(sy2inter,y2normreshape);

%% Rename to save typing later on

%Original data
x1 = series1x; %original data series 1x - not rounded
x2 = x1reshape; %factor interpolated data series 1x
x3 = series2x; %original data series 1x - not rounded
x4 = x2reshape; %factor interpolated data series 2x

y01 = series1y; %original series 1
y02 = series2y; %original series 2
y01n = series1ynorm; %standardised series 1
y02n = series2ynorm; %standardised series 2
y1ni = s1norminter; %standardised + interpolated series 2
y2ni = s2norminter; %standardised + interpolated series 2


%% PART 3 & 4: Plots, periodograms & peak identication 

% Plot the normalised and interpolated data 
% plot series 1 data with linear (red line plot) and then pchip (blue line plot) interpolations
% with datapoints plotted as black (k) filled circles (o)

figure(31)
ax(1) = subplot(2,1,1);
plot(x1,y01n,'.','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',2,'LineWidth',0.2), hold on
plot(x2,y1ni, 'b-','LineWidth',1), hold off
xlabel(xlbl1{1});
ylabel(ylbl1{1});
title(['\fontsize{12}',num2str(meanyr1x),'-yr PCHIP Interpolated Z-scores ' seriesname1{1}]);
xlim auto
ylim auto %([min(-6) max(4)])
%ylim auto %([min(-1) max(1)]);

%plot series 2 data with linear and pchip interpolations
ax(2) = subplot(2,1,2);
plot(x3,y02n,'k.','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',2,'LineWidth',0.2), hold on
plot(x4,y2ni, 'r-','LineWidth',1), hold off
xlabel(xlbl2{1});
ylabel(ylbl2{1});
title(['\fontsize{12}',num2str(meanyr2x),'-yr PCHIP Interpolated Z-scores ' seriesname2{1}]);
xlim auto
ylim auto %([min(-6) max(4)])
%ylim auto %([min(-1) max(1)]);

%link parts and zoom in to show the changes
linkaxes(ax(1:2),'x');
%xlim([min(series1x) max(series1x)])
xlim auto
ylim auto

h31 = figure(31)
set(h31, 'Units', 'centimeters','PaperSize', [20,20]);
print(h31, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig31', '-dpdf')

%% PART 3: Detrend the y data 

% y = detrend(series,'constant')
% constant detrends the linear trend in the time series by subtracting data 
% mean value from vector y series data and normalised y series data 

format long %remove detrend code here if dont want to detrend data first
series1yd =  detrend(series1y,'constant'); %series1y;
series2yd =  detrend(series2y,'constant'); %series2y;

series1ydnorm =  detrend(series1ynorm,'constant'); %series1ynorm;
series2ydnorm =  detrend(series2ynorm,'constant'); %series2ynorm;

% y = detrend(series,'linear')
% constant detrends the linear trend in the time series by subtracting data 
% mean value from vector y series data
format long
series1ydl =  detrend(series1y,'linear'); %series1y;
series2ydl =  detrend(series2y,'linear'); %series2y;

series1ydlnorm =  detrend(series1ynorm,'linear'); %series1ynorm;
series2ydlnorm =  detrend(series2ynorm,'linear'); %series2ynorm;

%% Write detrended data to single column text file
%fileID1 = fopen('Age (cal yr BP)', 'Series 1 - detrended');
%fprintf(fileID1, series1x, series1yd);
%fileID2 = fopen('Age (cal yr BP)', 'Series 2 - detrended');
%fprintf(fileID2, series1x, series1yd);

detrended1y = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/series1ydl.csv';
xlswrite(detrended1y,series1ydl);

detrended2y = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/series2ydl.csv';
xlswrite(detrended2y,series2ydl);

detrended1ynorm = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/series1ydlnorm.csv';
xlswrite(detrended1ynorm,series1ydlnorm);

detrended2ynorm = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/series2ydlnorm.csv';
xlswrite(detrended2ynorm,series2ydlnorm);

%% Plot detrended data 
% series1yd %this is long but can copy from data window if needed elsewhere

figure(32)
ax(1) = subplot(2,1,1);
plot(series1x,series1y,series1x,series1yd,series1x,series1ydl)
hold on
%plot(x1,y1,'.'); %check these are the same by running this as well
title(['\fontsize{12}','Detrended time series ' seriesname1{1}]);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
legend('Measured','Detrend (mean)','Detrend (LBF)');
legend('Location','eastoutside')
legend ('boxon');

ax(2) = subplot(2,1,2);
plot(series2x,series2y,series2x,series2yd,series2x,series2ydl)
hold on
%plot(x3,y3,'.'); %check these are the same by running this as well
title(['\fontsize{12}','Detrended time series ' seriesname2{1}]);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
legend('Measured','Detrend (mean)','Detrend (LBF)');
legend('Location','eastoutside')
legend ('boxon');

%link parts and zoom in to show the changes
linkaxes(ax(1:2),'x');
%xlim([min(series1x) max(series1x)])
xlim auto
ylim auto

h32 = figure(32)
set(h32, 'Units', 'centimeters','PaperSize', [20,20]);
print(h32, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig32', '-dpdf')

%% Plot detrended Z-scores data 

figure(33)
ax(1) = subplot(2,1,1);
plot(series1x,series1ynorm,series1x,series1ydnorm,series1x,series1ydlnorm)
hold on
%plot(x1,y1,'.'); %check these are the same by running this as well
title(['\fontsize{12}','Detrended time series Z-scores ' seriesname1{1}]);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
legend('Z-scores','Z-scores Dmean','Z-scores DLBF');
legend('Location','eastoutside')
legend ('boxon');

ax(2) = subplot(2,1,2);
plot(series2x,series2ynorm,series2x,series2ydnorm,series2x,series2ydlnorm)
hold on
%plot(x1,y1,'.'); %check these are the same by running this as well
title(['\fontsize{12}','Detrended time series Z-scores ' seriesname2{1}]);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
legend('Z-scores','Z-scores Dmean','Z-scores DLBF');
legend('Location','eastoutside')
legend ('boxon');

%link parts and zoom in to show the changes
linkaxes(ax(1:2),'x');
%xlim([min(series1x) max(series1x)])
xlim auto
ylim auto

h33 = figure(33)
set(h33, 'Units', 'centimeters','PaperSize', [20,20]);
print(h33, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig33', '-dpdf')

%% Interpolate detrended data using linear & PCHIP

%PCHIP is piecewise Cubic Hermite Interpolating Polynomial
%PCHIP is a shape preserving piecewise cubic interpolation which avoids
%artifects of splines and preserves shape of original data series
%better than linear interpolation

%Linear - dont use 
%series1L = interp1(series1x,series1yd,t1,'linear');
%series2L = interp1(series2x,series2yd,t2,'linear');

%PCHIP - mean
series1p = interp1(series1x,series1yd,t1,'pchip'); %t1 for adjusted time series = x2
series2p = interp1(series2x,series2yd,t2,'pchip'); %t1 for adjusted time series = x2

% PCHIP - linear - factor-yr datatset
series1pl = interp1(series1x,series1ydl,t1,'pchip'); %t2 for adjusted time series = x4
series2pl = interp1(series2x,series2ydl,t2,'pchip'); %t2 for adjusted time series on x4
series1plnorm = interp1(series1x,series1ydlnorm,t1,'pchip'); %t2 for adjusted time series = x4
series2plnorm = interp1(series2x,series2ydlnorm,t2,'pchip'); %t2 for adjusted time series on x4

% PCHIP - linear - 1-year dataset
%series1pl1 = interp1(series1x,series1yd,t3,'pchip'); %t2 for adjusted time series = x4
series2pl1 = interp1(series2x,series2yd,t4,'pchip'); %t2 for adjusted time series on x4
series1pl1norm = interp1(series1x,series1ydlnorm,t3,'pchip'); %t2 for adjusted time series = x4
series2pl1norm = interp1(series2x,series2ydlnorm,t4,'pchip'); %t2 for adjusted time series on x4

%% List of renamed variables to use as shortened versions from this point onwards

%Original data
x1 = series1x; %original data series 1x - not rounded
x2 = t1; %factor interpolated data series 1x
x3 = series2x; %original data series 1x - not rounded
x4 = t2; %factor interpolated data series 2x

y01 = series1y; %original data series 1
y02 = series2y; %original data series 2

%Detrended data
y1 = series1yd; %detrended data series 1 - subtracted from mean
y1l = series1ydl; %detrended data series 1 - subtracted from linear best fit
y1ln = series1ydlnorm; %normalised, detrended data series 1 - subtracted from linear best fit

y3 = series2yd; %detrended data series 2 - subtracted from mean
y3l = series2ydl; %detrended data series 2 - subtracted from linear best fit
y3ln = series2ydlnorm; %normalised,detrended data series 2 - subtracted from linear best fit

%Detrended & interpolated data
y2 = series1p; %detrended and factor interpolated data series 1 - subtracted from mean
y2l = series1pl; %detrended data series 1 - subtracted from linear best fit
y4 = series2p; %detrended and factor interpolated data series 2 - subtracted from mean
y4l = series2pl; %detrended and factor interpolated data series 2 - subtracted from linear

%Standardized, Detrended & interpolated data
y2ln = series1plnorm; %std and detrended data series 1 - subtracted from linear best fit
y4ln = series2plnorm; %stdized and detrended and factor interpolated data series 2 - subtracted from linear

%reshaped interpolated detrended y series into a vector with n rows x 1
%column from n columns x 1 row - to use later on for regression etc.

rs = 715

%Series 1
x2s = reshape(x2,[rs,1]);
y2s = reshape(y2,[rs,1]);
y2sl = reshape(y2l,[rs,1]);
y2sln = reshape(y2ln,[rs,1]);

%Series 2
x4s = reshape(x4,[rs,1]);
y4s = reshape(y4,[rs,1]);
y4sl = reshape(y4l,[rs,1]);
y4sln = reshape(y4ln,[rs,1]);

%Nomalise outputs - if needed
%y2slnorm = (y2sl-mean(y2sl))/std(y2sl);
%y4slnorm = (y4sl-mean(y4sl))/std(y4sl);

%write these values to csv - example - or copy paste from command window 
%need to open a blank excel spreadsheet first for this to work

%series 1
interx2s = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/x2s.csv';
xlswrite(interx2s,x2s);

intery2s = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/y2s.csv';
xlswrite(intery2s,y2s);
intery2sl = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/y2sl.csv';
xlswrite(intery2sl,y2sl);
intery2sln = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/y2sln.csv';
xlswrite(intery2sln,y2sln);
%interpolatedy2slnorm = 'y2slnorm.csv';
%xlswrite(interpolatedy2slnorm,y2slnorm);

%series 2
interx4s = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/x4s.csv';
xlswrite(interx4s,x4s);

intery4s = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/y4s.csv';
xlswrite(intery4s,y4s);
intery4sl = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/y4sl.csv';
xlswrite(intery4sl,y4sl);
intery4sln = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/y4sln.csv';
xlswrite(intery4sln,y4sln);
%interpolatedy4slnorm = 'y4slnorm.csv';
%xlswrite(interpolatedy4slnorm,y4slnorm);

%% Plot the detrended and interpolated data 

% plot series 1 data with linear (red line plot) and then pchip (blue line plot) interpolations
% with datapoints plotted as black (k) filled small circles (. or o) orange = #D95319,
% light blue = #4DBEEE, dk red  = #8b0000, darkgrey = #a9a9a9, darkblue = #00008b

figure (34)
ax(1) = subplot(4,1,1);
plot(x1,y1l,'.','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',2,'LineWidth',0.2), hold on
plot(x2,y2l, '-','Color','#00008b','LineWidth',1), hold off
%plot(x2,y2lst, 'r-','LineWidth',1), hold on
xlabel(xlbl1{1});
ylabel(ylbl1{1});
title(['\fontsize{12}',num2str(meanyr1x),'-yr PCHIP Interpolated/Detrended (LBF) ' seriesname1{1}]);
xlim auto
ylim auto %([min(-1) max(1)]);

%plot series 2 data with linear and pchip interpolations
ax(2) = subplot(4,1,2);
plot(x3,y3l,'.','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',2,'LineWidth',0.2), hold on
plot(x4,y4l, '-','Color','#8b0000','LineWidth',1), hold off
xlabel(xlbl2{1});
ylabel(ylbl2{1});
title(['\fontsize{12}',num2str(meanyr2x),'-yr PCHIP Interpolated/Detrended (LBF) ' seriesname2{1}]);
xlim auto
ylim auto %([min(-1) max(1)]);

%plot series 1 and 2 standardized data with linear pchip interpolation
ax(3) = subplot(4,1,3);
plot(x1,y1ln,'.','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',2,'LineWidth',0.2), hold on
plot(x2,y2ln, '-','Color','#00008b','LineWidth',1), hold off
xlabel(xlbl1{1});
ylabel(ylbl1{1});
title(['\fontsize{12}',num2str(meanyr1x),'-yr PCHIP Interpolated/Detrended Z-scores (LBF) ' seriesname1{1}]);
xlim auto
ylim auto %([min(-1) max(1)]);

ax(4) = subplot(4,1,4);
plot(x3,y3ln,'k.','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',2,'LineWidth',0.2), hold on
plot(x4,y4ln, '-','Color','#8b0000','LineWidth',1), hold off
xlabel(xlbl2{1});
ylabel(ylbl2{1});
title(['\fontsize{12}',num2str(meanyr2x),'-yr PCHIP Interpolated/Detrended Z-scores (LBT) ' seriesname2{1}]);
xlim auto
ylim auto %([min(-1) max(1)]);

h34 = figure(34)
set(h34, 'Units', 'centimeters','PaperSize', [20,20]);
print(h34, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig34', '-dpdf')

%% Scatterplot of interpolated-year detrended data compared with 10 or 100 year interpolated data 
figure (35)
ax(1) = subplot(2,1,1)
%scatter(series1pl1,series2pl1, '.','MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',0.2);
hold on;
scatter(y2l,y4l, 'o','MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',0.2);xlabel(ylbl1{1});
hold off;
ylabel(ylbl2{1});
title(['\fontsize{12}', seriesname1{1}, ' vs ',seriesname2{1}, ' PCHIP detrended']);
legend('10-yr Z-DLBF');
legend('Location','eastoutside')
legend ('boxoff');
xlim auto
ylim auto

ax(2) = subplot(2,1,2)
%scatter(series1plnorm,series2plnorm, '.','MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',0.2);
hold on;
scatter(y2ln,y4ln, 'o','MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',0.2);xlabel(ylbl1{1});
hold off;
ylabel(ylbl2{1});
title(['\fontsize{12}', seriesname1{1}, ' vs ',seriesname2{1},' PCHIP,Std & detrended']);
legend('10-yr Z-DLBF');
legend('Location','eastoutside')
legend ('boxoff');
xlim auto
ylim auto

h35 = figure(35)
set(h35, 'Units', 'centimeters','PaperSize', [20,20]);
print(h35, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig35', '-dpdf')

%% Frequency and periodgram plots for pchip data using Periodogram function

%Graphs in this section make the function for the frequency power spectrum 
%using the from pchip and detrended and interpolated data as input (y2ln
%for series 1 and y4ln for series 2
%256 datapoints and frequency of 1/mean interval value (interpolated 10
%yr data)
%512 and 1024 etc will produce smoother and more high density of peaks
%need to experiment to get the best combination for your data

%then plot 1/f as this is equal to the periodicity in years in the time series 
%uses rounded 10-yr interpolated data as input
%factor = 10
%check that meanyr1x and meanyr2 = 10

%series 1 - interpolated frequency plot
figure(36)
ax(1) = subplot(2,2,1);
[Pxx1,f1] = periodogram(y2ln,[],256,1/meanyr1x);
plot(f1,Pxx1,'k');
xlim ([min(0) max(0.04)]) %check whole range using 'xlim auto' first and then rescale
ylim auto %([min(0) max(200)]) %check whole range using 'ylim auto' first and then rescale
xlabel('Frequency')
ylabel('Power')
title(['\fontsize{10}','Auto-spectrum of ' seriesname1{1}]);

%series 1 - interpolated years plot 
ax(2) = subplot(2,2,2);
[Pxx1,f1] = periodogram(y2ln,[],256,1/meanyr1x);
semilogx(1./f1,Pxx1,'k');
xlim ([min(0) max(10000)])
xticks auto %use auto first to check range
ylim auto %([min(0) max(200)]) %use auto first to check range
xlabel('Years')
ylabel('Power')
title(['\fontsize{10}','Periodicity of ' seriesname1{1}]);
hold on

% peak identification to find all the data points as the peaks
MPPf1= std(Pxx1)/4;
%stdev of Pxx1 to pick out major peaks >2*stdev (2s, 95.4%) of power
%spectrum
[pksf1,locsf1] = findpeaks(Pxx1,'MinPeakProminence',MPPf1);
% offset values of peak heights for plotting
plot(1./f1(locsf1),pksf1+5,'v','markerfacecolor','b','markeredgecolor','b', 'MarkerSize',4);
hold off;
%location of series 1 signifcant peaks in years and peak value
period1 = 1./f1(locsf1)
pksf1

%link y axes for these plots 
%linkaxes([ax(1),ax(2)],'y')

%series 2 - interpolated frequency plot 
ax(3) = subplot(2,2,3);
[Pxx2,f2] = periodogram(y4ln,[],248,1/meanyr2x);

plot(f2,Pxx2,'k');
xlim ([min(0) max(0.04)]) %check whole range using 'xlim auto' first
ylim auto %([min(0) max(200)]) %check whole range using 'ylim auto' first
xlabel('Frequency')
ylabel('Power')
title(['\fontsize{10}','Auto-spectrum of ' seriesname2{1}]);

%series 2 - interpolated years plot
ax(4) = subplot(2,2,4);
[Pxx2,f2] = periodogram(y4ln,[],248,1/meanyr2x);
semilogx(1./f2,abs(Pxx2),'k')
xlim ([min(0) max(10000)]) %use auto first to check range
xticks auto
ylim auto %([min(0) max(200)]) %use auto first to check range
xlabel('Years')
ylabel('Power')
title(['\fontsize{10}','Periodicity of ' seriesname2{1}]);
hold on

% use peak identification to find all the data points as the peaks
MPPf2 = std(abs(Pxx2))/4;
%stdev of Pxx1 to pick out major peaks >2*stdev (2s, 95.4%) of power
%spectrum
[pksf2,locsf2] = findpeaks(Pxx2,'MinPeakProminence',MPPf2);
% offset values of peak heights for plotting
plot(1./f2(locsf2),pksf2+5,'v','markerfacecolor','b','markeredgecolor','b', 'MarkerSize',4);
hold off;
%location of series 2 signifcant peaks in years and peak value
period2 = 1./f2(locsf2)
pksf2

%use linkaxes to link y axes for series 2 plots - if they have similar axes values 
%linkaxes([ax(3),ax(4)],'y')

%link axes for series 1 and series 2 x-axis (log years) plots
%linkaxes([ax(2),ax(4)],'x')

%Write peak locations to file - s1 = series 1 periodicity cycle length in
%years -then compare values to s2 to see if there are simialrities - group
%<100 year cycles into eg 30-50 years if there are lots of decadal scale cycles 
s1period = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/s1period.csv';
xlswrite(s1period,period1);

s2period = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/s2period.csv';
xlswrite(s2period,period2);

h36 = figure(36)
set(h36, 'Units', 'centimeters','PaperSize', [20,20]);
print(h36, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig36', '-dpdf')

%% Magnitude squared coherence & Frequency Cross Spectrum of the time series (MSC)
% This plot identifies frequency domain correlations between two time series.
% Phase estimates in the cross-spectrum are where significant frequency-
% domain correlation exists - mcohere (hamming type) uses Welch's overlapped
% segment averaging (WOSA) and multitaper techniques with a window length
% of 100 samples, equivalent to 10 periods of a 100 Hz signal and overlap
% of 80 - therefore need to mulitply the meanyr interpolated interval by
% interpolation factor to give an output that is in Hz

figure (37)
ax(1) = subplot(2,2,1);
[Cxy,f4]= mscohere(y2ln,y4ln,hamming(100),80,meanyr1x*(factor1)^2);
plot(f4,Cxy,'k');
xlabel('Frequency (Hz)');
ylabel('Magnitude Squared Coherence');
title(['\fontsize{8}','Coherence of ' seriesname1{1}, ' & ',seriesname2{1}]);
xlim auto %([min(0) max(0.1)])
ylim auto %([min(0) max(1)])
f4

ax(2) = subplot(2,2,2);
[Cxy,f4]= mscohere(y2ln,y4ln,hamming(100),80,meanyr1x*(factor2)^2);
% x-axis in years
xf4 = (1./f4)*factor1^2;
semilogx(xf4,Cxy,'k');
xlabel('Years')
ylabel('Magnitude Squared Coherence')
title(['\fontsize{8}','Coherence of ' seriesname1{1}, ' and ',seriesname2{1}]);
xlim auto %([min(0) max(10000)]);
xticks auto %([0:200:1000]);
ylim ([min(0) max(1)])
hold on

% peak identification to find signifcant series 1 peaks
% defined as promiment peaks >2s offset in y axis value
MPPf4 = std(Cxy)/4;
%stdev of Pxx1 to pick out major peaks >2*stdev (2s, 95.4%) of power
%spectrum
[pksf4,locsf4] = findpeaks(Cxy,'MinPeakProminence',MPPf4);
% offset values of peak heights for plotting
plot(xf4(locsf4),pksf4+0.1,'v','markerfacecolor','b','markeredgecolor','b', 'MarkerSize',4);
hold off;
%location of coherent peaks in years and peak value
mscpeak = xf4(locsf4)
pksf4
%link y axes for Series 1 plots
linkaxes([ax(1),ax(2)],'y')

%Write msc peak locations to file 
if isempty(pksf4)
    disp('No peak locations')
else
    mscpk = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/msc.csv';
    xlswrite(mscpk,mscpeak);
end

% Frequency cross spectrum 
ax(3) = subplot(2,2,3);
[Pxy,f3]= cpsd(series1pl,series2pl,hamming(100),80,meanyr1x*factor1);
plot(f3,abs(Pxy),'k');
xlabel('Frequency (Hz)')
ylabel('Power')
xlim auto %([min(0) max(0.1)])
ylim auto %([min(0) max(1.2)])
title(['\fontsize{8}','Cross-Spectrum of ' seriesname1{1}, ' & ',seriesname2{1}]);
hold on

ax(4) = subplot(2,2,4);
[Pxy,f3]= cpsd(series1pl,series2pl,hamming(100),80,meanyr1x*factor1);
%create new axis from factor^2 = 100
xf3 = (1./f3)*factor1^2; 
semilogx(xf3,abs(Pxy),'k');
xlabel('Years')
ylabel('Power')
title(['\fontsize{8}','Cross-Spectrum of ' seriesname1{1}, ' & ',seriesname2{1}]);
xlim auto %([min(0) max(10000)]);
xticks auto %([0:100:1000]);
ylim auto;
hold on

% peak identification to find significant Series 2 peaks
% defined as promiment peaks >2s offset in y axis
MPPf3 = std(abs(Pxy))/4;
%stdev of Pxx1 to pick out major peaks >2*stdev (2s, 95.4%) of power
%spectrum
[pksf3,locsf3] = findpeaks(abs(Pxy),'MinPeakProminence',MPPf3);
% offset values of peak heights for plotting
plot(xf3(locsf3),pksf3+0.01,'v','markerfacecolor','b','markeredgecolor','b', 'MarkerSize',4);
hold off;
%location of cross spectrum peaks in years and peak value
xf3(locsf3)
pksf3

%Write frcs peak locations to file 
if isempty(pksf3)
    disp('No peak locations')
else
    frcspect = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/frcspect.csv';
    xlswrite(frcspect,xf3(locsf3));
end

%link y axes for series 2 plots 
%linkaxes([ax(3),ax(4)],'y')
%link axes for series 1 and series 2 x-axis (log years) plots
linkaxes([ax(2),ax(4)],'x')

%alternate way of resizing x axis - but cant list years as output 
% xticks([0:2:10]) 
% title(['\fontsize{8}','Cross-Spectrum of ' seriesname1{1}, ' & ',seriesname2{1}]);
% Resize x axis to plot in years need to x100 
% xt = get(gca, 'XTick'); 
% Relabel 'XTick' With 'XTickLabel' Values
% set(gca, 'XTick', xt, 'XTickLabel', (xt*factor^2))
% needs to by x100 when each series has a factor = 10, i.e., factor^2

%Save 4-part figure to file
h37 = figure(37)
set(h37, 'Units', 'centimeters','PaperSize', [20,20]);
print(h37, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig37', '-dpdf')

%% Phase lag (in radians) of cross spectrum
% plot is of the phase of the cross spectrum
% which indicates frequencies with sig coherence between two series
% and shows known phase between their sinusoidal components of input
% signals

figure(38)
phase = angle(Pxy);
plot(f3,-angle(Pxy)/pi);
xlabel('Frequency (Hz)')
ylabel('Lag \times\pi rad')
title(['\fontsize{12}','Cross spectrum Phase of ' seriesname1{1}, ' & ',seriesname2{1}]);
xlim auto %([min(0) max(0.1)])
ylim auto %([min(-4) max(4)])

h38 = figure(38)
set(h38, 'Units', 'centimeters','PaperSize', [20,20]);
print(h38, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig38', '-dpdf')

%% Evolutionary power spectrum - spectrogram map
%Computes a short-time Forier transform (STFT) of overlapping segments of the
%time series using a Hamming window of 400 datapoints, 50 data point overlap 
% STFT 256 resolution (higher produces a smoother plot) and frequency set
% to sample as above 
% not sure if this is correct - need better interpretation of what this
% means - check axes are labelled properly 

figure(39)
spectrogram(series1pl,64,50,256,1/factor1);
title(['\fontsize{12}','Evolutionary Power Spectrum of ' seriesname1{1}]);
xlabel('Frequency (1/yr)')
ylabel('Time (yr)')
colormap jet;
xlim auto %([min(0) max(150)]);
ylim auto %([min max]);

h39 = figure(39)
set(h39, 'Units', 'centimeters','PaperSize', [20,20]);
print(h39, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig39', '-dpdf')

figure(40)
spectrogram(series2pl,64,50,256,1/factor1);
title(['\fontsize{12}','Evolutionary Power Spectrum of ' seriesname2{1}]);
xlabel('Frequency (1/yr)')
ylabel('Time (yr)')
colormap jet;
xlim auto %([min(0) max(150)]);
ylim auto %([min max]);

h40 = figure(40)
set(h40, 'Units', 'centimeters','PaperSize', [20,20]);
print(h40, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig40', '-dpdf')

%% Least squares polynomial fitting
% Plot series 1 data and compare to interpolated plot
% dk red  = #8b0000, darkgrey = #a9a9a9, darkblue = #00008b

figure(41)
ax(1) = subplot(4,1,1);
plot(x1,series1ydl,'.','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',2,'LineWidth',0.2), hold on
plot(x2,y2l,'-','Color', '#00008b', 'LineWidth',1), hold off
xlabel(xlbl1{1});
ylabel(ylbl1{1});
title(['\fontsize{12}',num2str(meanyr1x),'-yr PCHIP Interpolated/Detrended (LBT) ' seriesname1{1}]);
axis auto %([yrmin1 yrmax1 min(y1) max(y1)]);
ylim auto %([min(-1) max(1)]);
xlim auto %([min(x1) max(x1)]);
%yticks ([-10 -5 -2 -1 0 1 2 5 10])
 
ax(2) = subplot(4,1,2);
plot(x3,series2ydl,'.','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',2,'LineWidth',0.2), hold on
plot(x4,y4l,'-', 'Color', '#8b0000', 'LineWidth',1), hold off
xlabel(xlbl2{1});
ylabel(ylbl2{1});
title(['\fontsize{12}',num2str(meanyr2x),'-yr PCHIP Interpolated/Detrended (LBT) ' seriesname2{1}]);
axis([yrmin1 yrmax1 min(y1) max(y1)]);
xlim auto
ylim auto %([min(-1) max(1)]);
xlim auto %([min(x1) max(x1)]);
%yticks ([-10 -5 -2 -1 0 1 2 5 10])

ax(3) = subplot(4,1,3);
plot(x1,series1ydlnorm,'.','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',2,'LineWidth',0.2), hold on
plot(x2,y2ln,'-','Color', '#00008b', 'LineWidth',1), hold off
xlabel(xlbl1{1});
ylabel(ylbl1{1});
title(['\fontsize{12}',num2str(meanyr1x),'-yr PCHIP Interpolated/Detrended (LBT) Z-scores ' seriesname1{1}]);
axis auto %([yrmin1 yrmax1 min(y1) max(y1)]);
ylim auto %([min(-1) max(1)]);
xlim auto %([min(x1) max(x1)]);
yticks ([-10 -5 -2 -1 0 1 2 5 10]) %auto %([0:100:1000]);
 
ax(4) = subplot(4,1,4);
plot(x3,series2ydlnorm,'.','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',2,'LineWidth',0.2), hold on
plot(x4,y4ln,'-', 'Color', '#8b0000', 'LineWidth',1), hold off
xlabel(xlbl2{1});
ylabel(ylbl2{1});
title(['\fontsize{12}',num2str(meanyr2x),'-yr PCHIP Interpolated/Detrended (LBT) Z-scores ' seriesname2{1}]);
axis([yrmin1 yrmax1 min(y1) max(y1)]);
xlim auto
ylim auto %([min(-1) max(1)]);
xlim auto %([min(x1) max(x1)]);
yticks ([-10 -5 -2 -1 0 1 2 5 10])

linkaxes([ax(1),ax(2), ax(3), ax(4)], 'x');
linkaxes([ax(1),ax(2), ax(3), ax(4)], 'y');

h41 = figure(41)
set(h41, 'Units', 'centimeters','PaperSize', [20,20]);
print(h41, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig41', '-dpdf')

%% Least squares polynomial fitting
% four part figure showing original data and residuals in plot 1 and 2
% and detrended and interpolated data and its residuals in plot 3 and 4
% this shows the effect of removing the linear trend from the original
% time series with DW increasing if it has been removed

% The Durbin-Watson test assesses whether there is autocorrelation
% among the residuals or not.

% For example - value of the Durbin-Watson test statistic is 2.0526. 
% The p-value of 0.6285 suggests that the residuals are not autocorrelated.
% p = dwtest(r,x) returns the p-value for the Durbin-Watson test of the 
% null hypothesis that the residuals from a linear regression are uncorrelated. 
% The alternative hypothesis is that there is autocorrelation among the residuals.
% A significantly small p-value casts doubt on the validity of the null hypothesis 
% and indicates correlation among residuals.

%Series 1 plots 

figure(42)
ax(1) = subplot(4,1,1)
plot(x1,y01)
hold on
title(['\fontsize{10}Least squares polynomial fit for ', seriesname1{1}]);
%linear regression model
mdl1 = fitlm(x1,series1y);
r1sq = mdl1.Rsquared.adjusted
mdl1

% least squares polynomial fit
% p = polyfit(x,y,n) returns the coefficients for a polynomial p(x) of 
% degree n that is a best fit (in a least-squares sense) for the data in y
coeffs = polyfit(x1,series1y,1);
%fits least squares polynomial linear regression line to data
ylabel(ylbl1{1});
%fit linear curve to interpolated data
vq1fit = coeffs(2)+coeffs(1)*x1; 
plot(x1,vq1fit,'linewidth',1)
residuals1 = series1y - vq1fit;

title(['\fontsize{10}Least squares polynomial fit: ', seriesname1{1}]);
xlim auto
ylim auto

ax(2) = subplot(4,1,2)
[b1,bint1,r1] = regress(y01,x1);
plot(x1,residuals1,'linewidth',1);%residuals from fit above
ylabel('Residuals');
xlim auto
ylim auto
format long
[p1,DW1] = dwtest(r1,x1) %exact p values in two tailed test
title(['\fontsize{10}Residuals p=',num2str(p1),'; DW=',num2str(DW1)])

ax(3) = subplot(4,1,3)
plot (x2s,y2sl);
hold on
coeffs = polyfit(x2,y2,1);%fits least squares polynomial linear regression line to data
ylabel(ylbl1{1});
%fit linear curve to interpolated data
vq2fit = coeffs(2)+coeffs(1)*x2;
title(['\fontsize{10}Least squares polynomial fit: ',num2str(meanyr1x),'-yr PCHIP Detrended ' seriesname1{1}]); 
plot(x2,vq2fit,'linewidth',1), hold off;
xlim auto
ylim auto

%reshape the interpolated detrended y series into a vector with n rows x 1
%column from n columns x 1 row - to use later on for regression etc.
x2s = reshape(x2,[rs,1])
y2s = reshape(y2,[rs,1])
x4s = reshape(x4,[rs,1])
y4s = reshape(y4,[rs,1])

x2s = reshape(x2,[rs,1])
y2sl = reshape(y2l,[rs,1])
x4s = reshape(x4,[rs,1])
y4sl = reshape(y4l,[rs,1])

% regress and plot the data
ax(4) = subplot(4,1,4);
residuals2 = y2sl - vq2fit;
[b2,bint2,r2] = regress(y2s,x2s);
plot(x2s,residuals2,'linewidth',1)%residuals from fit above
xlabel(xlbl1{1});
ylabel('Residuals');
xlim auto
ylim auto
format long
[p2,DW2] = dwtest(r2,x2s) %exact p values in two tailed test
title(['\fontsize{10}Residuals p=',num2str(p2),'; DW=',num2str(DW2)])

% Write output summary stats to file
diary('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/mdl1.txt')
mdl1
p1
DW1
p2
DW2
diary off

h42 = figure(42)
set(h42, 'Units', 'centimeters','PaperSize', [20,20]);
print(h42, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig42', '-dpdf')

%% Summary stats output for Series 1
% write these values to file
% filenamestats1 = 'Polyfit_stats_Series1.csv';
% series1stats = {'p-value','DW stat'; p2, DW2}; 
% xlswrite(filenamestats1, series1stats);

%% Polyfit Series 2
%Series 2 plots 

figure(43)
ax(1) = subplot(4,1,1);
plot(series2x,series2y)
hold on

%linear regression model
mdl2 = fitlm(x3,y02);
r2sq = mdl2.Rsquared.adjusted
mdl2

% least squares polynomial fit
% p = polyfit(x,y,n) returns the coefficients for a polynomial p(x) of 
% degree n that is a best fit (in a least-squares sense) for the data in y
coeffs = polyfit(x3,y02,1);
%fits least squares polynomial linear regression line to data
ylabel(ylbl2{1});
%fit linear curve to interpolated data
vq3fit = coeffs(2)+coeffs(1)*series2x; 
plot(x3,vq3fit,'linewidth',1)
residuals3 = y02 - vq3fit;
title(['\fontsize{10}Least squares polynomial fit for ', seriesname2{1}]);
xlim auto
ylim auto %([min(series2y) max(series2y)]) %auto

ax(2) = subplot(4,1,2);
[b3,bint3,r3] = regress(y02,x3);
plot(x3,residuals3,'linewidth',1);%residuals from fit above
ylabel('Residuals');
xlim auto
ylim auto
format long
[p3,DW3] = dwtest(r3,x3) %exact p values in two tailed test
title(['\fontsize{10}Residuals p=',num2str(p3),'; DW=',num2str(DW3)])

ax(3) = subplot(4,1,3);
plot (x4s,y4sl);
hold on
coeffs = polyfit(x4s,y4sl,1);%fits least squares polynomial linear regression line to data
ylabel(ylbl2{1});
%fit linear curve to interpolated data
vq4fit = coeffs(2)+coeffs(1)*x4s;
title(['\fontsize{10}Least squares polynomial fit for ',num2str(meanyr2x),'-yr PCHIP Interpolation of Detrended ' seriesname2{1}]); 
plot(x4s,vq4fit,'linewidth',1), hold off;
residuals4 = y4sl - vq4fit;
xlim auto
ylim ([min(y4) max(y4)]) %auto

ax(4) = subplot(4,1,4);
[b4,bint4,r4] = regress(y4sl,x4s);
plot(x4,residuals4,'linewidth',1)%residuals from fit above
xlabel(xlbl2{1});
ylabel('Residuals');
xlim auto
ylim auto
format long
[p4,DW4] = dwtest(r4,x4s) %exact p values in two tailed test
title(['\fontsize{10}Residuals p=',num2str(p4),'; DW=',num2str(DW4)])

linkaxes([ax(1),ax(2),ax(3),ax(4)],'x');

% Write output summary stats to file
diary('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/mdl2.txt')
mdl2
p3
DW3
p4
DW4
diary off

h43 = figure(43)
set(h43, 'Units', 'centimeters','PaperSize', [20,20]);
print(h43, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig43', '-dpdf')

%% Summary stats output for Series 2
% write these values to file
% filenamestats4 = 'Polyfit_stats_Series1.csv';
% series2stats = {'p-value','DW stat'; p4, DW4}; 
% csvwrite(filenamestats4, series2stats);

%% Serial Autocorrelation
% autocorrelation using residual plots
% using 'residuals' from polynomial trend fitting plot

figure(44)
%plot the autocorrelation for detrended and interpolated series 1
ax(1) = subplot(2,2,1);
acf1 = autocorr(y2l);
%x = autocorr(v1) for original data
semilogx(acf1);
xlabel('Lag time (years)');
ylabel('Correlation coefficient');
title(['\fontsize{10}AC of Detrend-Int ', seriesname1{1}]);

%plot the autocorrelation for derended and interpolated series 2
ax(2) = subplot(2,2,2);
acf2 = autocorr(y4l);
%x = autocorr(v1) for original data
semilogx(acf2);
xlabel('Lag time (years)');
ylabel('Correlation coefficient');
title(['\fontsize{10}AC of Detrend-Int ', seriesname2{1}]);

%plot the autocorrelation for detrended and interpolated series 1
ax(3) = subplot(2,2,3);
acf3 = autocorr(y2ln);
%x = autocorr(v1) for original data
semilogx(acf3);
xlabel('Lag time (years)');
ylabel('Correlation coefficient');
title(['\fontsize{10}AC of Detrend-Int Z-scores ', seriesname1{1}]);

%plot the autocorrelation for derended and interpolated series 2
ax(4) = subplot(2,2,4);
acf4 = autocorr(y4ln);
%x = autocorr(v1) for original data
semilogx(acf4);
xlabel('Lag time (years)');
ylabel('Correlation coefficient');
title(['\fontsize{10}AC of Detrend-Int Z-scores ', seriesname2{1}]);

h44 = figure(44)
set(h44, 'Units', 'centimeters','PaperSize', [20,20]);
print(h44, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig44', '-dpdf')

%% Peak identification

% Series 1
figure(45)

% A) peak identification to find all the data points as the peaks
ax(1) = subplot(4,1,1);
plot(x2,y2sln); hold on;
[pks0,locs0] = findpeaks(y2sln,x2,'MinPeakProminence',0);
plot(locs0,pks0+0.5,'bv','markerfacecolor',[0 0 1], 'MarkerSize',3);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
locs0yrs= locs0;
locs0mean = round(mean(diff(locs0yrs)))
locs0std = round(std(diff(locs0yrs)))
title(['\fontsize{10} All peaks (',num2str(meanyr1x),'-yr IDZ data): ',seriesname1{1},' Mean peak diff = ',num2str(locs0mean),'±',num2str(locs0std),' yrs'])

% B) 1s filtering - MinPeakProminence value = >stdev of interpolated data interval
MPP = std(intv1)
ax(2) = subplot(4,1,2);
plot(x2,y2sln); hold on;
[pks1,locs1] = findpeaks(y2sln,'MinPeakProminence',MPP);
% calculate offset values of peak heights for plotting
plot(x2(locs1),pks1+0.1,'bv','markerfacecolor',[0 0 1], 'MarkerSize',3);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
locs1yrs= x2(locs1)
locs1mean = round(mean(diff(locs1yrs)))
locs1std = round(std(diff(locs1yrs)))
title(['\fontsize{10} Prominent IDZ Peaks [Trough-Peak Diff >',num2str(MPP),' (1s)] Mean peak diff = ',num2str(locs1mean),'±',num2str(locs1std),' yrs'])%adds variable to title - remember []

% C) 2s filtering applied to all data
ax(3) = subplot(4,1,3);
[pks2,locs2] = findpeaks(y2sln,'MinPeakProminence',2*MPP);
plot(x2,y2sln); hold on;
% offset values of peak heights for plotting
plot(x2(locs2),pks2+0.1,'bv','markerfacecolor',[0 0 1], 'MarkerSize',3);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
locs2yrs= x2(locs2);
locs2mean = round(mean(diff(locs2yrs)))
locs2std = round(std(diff(locs2yrs)))
title(['\fontsize{10} IDZ peak-trough variance >2s Mean peak diff = ',num2str(locs2mean),'±',num2str(locs2std),' yrs']) %adds variable to title

% D) Saturated peaks MinPeakProminence value of >stdev of interpolated data interval/ T1 interpolation as a horizontal peak-width constraint
ax(4) = subplot(4,1,4);
plot(x2,y2sln); hold on;
TH = MPP/T1 %sets threshold value to stdev / interpolation years e.g., 10 or 50 years
[pks3,locs3] = findpeaks(y2sln,'threshold',TH);
% offset values of peak heights for plotting
plot(x2(locs3),pks3+0.1,'bv','markerfacecolor',[0 0 1], 'MarkerSize',3);
peakInterval3 = diff(locs3);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
locs3yrs= x2(locs3);
locs3mean = round(mean(diff(locs3yrs)))
locs3std = round(std(diff(locs3yrs)))
title(['\fontsize{10} IDZ w/o saturated peaks: Threshold > ',num2str(TH), ' (1s/int.) Mean peak diff = ',num2str(locs3mean),'±',num2str(locs3std),' yrs'])%adds variable to title

%link A-D and zoom in to show the changes
linkaxes(ax(1:4),'xy')
xlim([min(x1) max(x1)])
%ylim([0 max(v1-1)])
ylim auto

% Write output summary stats to file & save figure
diary('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Peaks_S1.txt')
locs0mean, locs0std, 
locs1mean, locs1std
locs2mean, locs2std
locs3mean, locs3std
locs0yrs, locs1yrs, locs2yrs, locs3yrs
diary off

h45 = figure(45)
set(h45, 'Units', 'centimeters','PaperSize', [20,20]);
print(h45, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig45', '-dpdf')

%% Series 2
figure(46)

% A) peak identification to find all the data points as the peaks
ax(1) = subplot(4,1,1);
plot(x4,y4sln); hold on;
[pks4,locs4] = findpeaks(y4sln,x4,'MinPeakProminence',0);
plot(locs4,pks4+0.5,'bv','markerfacecolor',[0 0 1], 'MarkerSize',3);
% offset values of peak heights for plotting
xlabel(xlbl2{1});
ylabel(ylbl2{1});
locs4yrs= locs4;
locs4mean = round(mean(diff(locs4yrs)))
locs4std = round(std(diff(locs4yrs)))
title(['\fontsize{10} All peaks (',num2str(meanyr2x),'-yr IDZ data)',seriesname2{1},' Mean peak diff = ',num2str(locs4mean),'±',num2str(locs4std),' yrs'])

% B) 1s filtering and caluclating distance between peaks in first plot
% uses MinPeakProminence value of >stdev of interpolated data interval
MPP = std(intv2) %stdev_all as alternative
ax(2) = subplot(4,1,2);
plot(x4,y4sln); hold on;
[pks5,locs5] = findpeaks(y4sln,'MinPeakProminence',MPP);
% calculate offset values of peak heights for plotting
plot(x4(locs5),pks5+0.1,'bv','markerfacecolor',[0 0 1], 'MarkerSize',3);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
locs5yrs= x4(locs5);
locs5mean = round(mean(diff(locs5yrs)))
locs5std = round(std(diff(locs5yrs)))
title(['\fontsize{10} Prominent IDZ Peaks [Trough-Peak Difference >',num2str(2*MPP),' (1s)] Mean peak diff = ',num2str(locs5mean),'±',num2str(locs5std),' yrs'])%adds variable to title - remember []
 
% C) 2s filtering applied to all data
ax(3) = subplot(4,1,3);
[pks6,locs6] = findpeaks(y4sln,'MinPeakProminence',2*MPP);
plot(x4,y4sln); hold on;
% offset values of peak heights for plotting
plot(x4(locs6),pks6+0.1,'bv','markerfacecolor',[0 0 1], 'MarkerSize',3);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
locs6yrs= x4(locs6);
locs6mean = round(mean(diff(locs6yrs)))
locs6std = round(std(diff(locs6yrs)))
title(['\fontsize{10} IDZ peak-trough variance >2s Mean peak diff = ',num2str(locs6mean),'±',num2str(locs6std),' yrs'])%adds variable to title

% D) Saturated peaks MinPeakProminence value of >stdev of interpolated data interval/ T1 interpolation as a horizontal peak-width constraint
ax(4) = subplot(4,1,4);
plot(x4,y4sln); hold on;
TH = (MPP)/T2 %sets threshold value to 2xstdev / interpolation value e.g., 10 or 50 years
[pks7,locs7] = findpeaks(y4sln,'threshold',TH);
% offset values of peak heights for plotting
plot(x4(locs7),pks7+0.1,'bv','markerfacecolor',[0 0 1], 'MarkerSize',3);
peakInterval4 = diff(locs7);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
locs7yrs= x4(locs7);
locs7mean = round(mean(diff(locs7yrs)))
locs7std = round(std(diff(locs7yrs)))
title(['\fontsize{10} IDZ w/o saturated peaks: Threshold > ',num2str(TH), ' (1s/int.) Mean peak diff = ',num2str(locs7mean),'±',num2str(locs7std),' yrs'])%adds variable to title
 
%link A-D and zoom in to show the changes
linkaxes(ax(1:4),'xy')
xlim([min(x2) max(x2)])
%ylim([0 max(v1-1)])
ylim auto

% Write output summary stats to file & save figure
diary('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Peaks_S2.txt')
locs4mean, locs4std, 
locs5mean, locs5std
locs6mean, locs6std
locs7mean, locs7std
locs4yrs, locs5yrs, locs6yrs, locs7yrs
diary off

h46 = figure(46)
set(h46, 'Units', 'centimeters','PaperSize', [20,20]);
print(h46, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig46', '-dpdf')
 
%% PEAK FILTERING

%Peak-prominence >mean+2s of detrended and interpolated data 
figure(47)
ax(1) = subplot(2,1,1);
%series 1
%measuring width of the peaks using half prominence as a reference
P1 = meanyr1+(2*stdevyr1); %sets MinPeakProminence value to >2stdev
findpeaks(y2,x2,'MinPeakDistance',P1); %'Annotate','extents','WidthReference','halfheight');
[pks8,locs8] = findpeaks(y2,x2,'MinPeakDistance',P1); %'Annotate','extents','WidthReference','halfheight');
hold on
%plot(y2,x2)
plot(locs8,pks8,'x')
xlabel(xlbl1{1});
ylabel(ylbl1{1});
xlim auto;
ylim auto;
locs8yrs = locs8*factor1;
locs8mean = mean(diff(locs8))
locs8std = std(diff(locs8))
title(['\fontsize{10}Peak-filtered ',num2str(meanyr1x), '-yr ID ',seriesname1{1},' >Mean+2s dist; Mean peak diff = ',num2str(round(locs8mean),-2),'±',num2str(round(locs8std),-2),' yrs']);

%series2
ax(2) = subplot(2,1,2);
P2 = meanyr2+(2*stdevyr2); %sets MinPeakProminence value to >2stdev
findpeaks(y4,x4,'MinPeakDistance',P2); %'Annotate','extents','WidthReference','halfheight');
[pks9,locs9] = findpeaks(y4,x4,'MinPeakDistance',P2); %'Annotate','extents','WidthReference','halfheight');
hold on
%plot(y2,x2)
plot(locs9,pks9,'x')
xlabel(xlbl1{1});
ylabel(ylbl1{1});
xlim auto;
ylim auto;
locs9yrs = locs9*factor1;
locs9mean = mean(diff(locs9))
locs9std = std(diff(locs9))
title(['\fontsize{10}Peak-filtered ',num2str(meanyr1x), '-yr ID ',seriesname2{1},' >Mean+2s dist; Mean peak diff = ',num2str(round(locs9mean),-2),'±',num2str(round(locs9std),-2),' yrs']);

% Write output summary stats to file & save figure
diary('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Peak_prom_2s.txt')
locs8mean, locs8std, 
locs9mean, locs9std
%locs8, locs9
diary off

h47 = figure(47)
set(h47, 'Units', 'centimeters','PaperSize', [20,20]);
print(h47, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig47', '-dpdf')

%% Peak identification using prominence >1/2 width 

figure(48)
%Peak prominence >1/2 width of standardised, detrended and interpolated peak prominence input data 
%series 1
[pks10,locs10] = findpeaks(y2sln,x2,'Annotate','extents');
findpeaks(y2sln,x2,'Annotate','extents');
legend('Location','southoutside');
xlabel(xlbl1{1});
ylabel(ylbl1{1});
xlim auto;
ylim auto;
locs10mean = mean(diff(locs10))
locs10std = std(diff(locs10))
title(['\fontsize{10}Peak-filtered ',num2str(meanyr1x), '-yr IDZ data ',seriesname1{1},' Peak width >1/2 height; Mean peak diff = ',num2str(round(locs10mean),-2),'±',num2str(round(locs10std),-2),' yrs']);

h48 = figure(48)
set(h48, 'Units', 'centimeters','PaperSize', [20,20]);
print(h48, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig48', '-dpdf')

figure(49)
%series
[pks11,locs11] = findpeaks(y4sln,x4,'Annotate','extents');
findpeaks(y4sln,x4,'Annotate','extents');
legend('Location','southoutside');
%P = stdev_all/2; %sets MinPeakProminence value to >stdev
%[a,b] = findpeaks(y2,x2,'MinPeakProminence',P)
%findpeaks(y1,x2,'MinPeakProminence',P)
xlabel(xlbl2{1});
ylabel(ylbl2{1});
xlim auto
ylim auto
locs11mean = mean(diff(locs11))
locs11std = std(diff(locs11))
title(['\fontsize{10}Peak-filtered ',num2str(meanyr2x), '-yr IDZ ' ,seriesname2{1},' Peak width >1/2 height; Mean peak diff = ',num2str(round(locs11mean),-2),'±',num2str(round(locs11std),-2),' yrs']);

%linkaxes([ax(1),ax(2)],'x')

% Write output summary stats to file & save figure
diary('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Peak_prom_0.5width.txt')
locs10mean, locs10std, 
locs11mean, locs11std
%locs10, locs11
diary off

h49 = figure(49)
set(h49, 'Units', 'centimeters','PaperSize', [20,20]);
print(h49, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig49', '-dpdf')

%% PART 5 - CHANGEPOINT ANALYSIS

% most significant changepoint - for interpolted dataset
% based on largest change in the root-mean-square
% after performing element-wise differentiation to remove any slowly 
% varying trends - then adds  zone added at end

% this runs on detrended series 1 and series 2 data overladi by interpolated factor 1 data
% to use standardized data replace
% y1l > y1ln - data points for series 1 detrended
% y2l > y2ln - data points for series 1 detrended & interpolated
% y3l > y1ln - data points for series 2 detrended
% y4l > y4ln - data points for series 2 detrended & interpolated 

figure(51)
% series 1
ax(1) = subplot(2,1,1);
plot(x1,y1ln,'ko','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',1,'LineWidth',0.2), hold on
plot(x2,y2ln,'b','LineWidth',1);
i = findchangepts(diff(y2),'Statistic','rms');
ax = gca;
xp1 = [x2(i) ax.XLim([2 2]) x2(i)];
yp1 = ax.YLim([1 1 2 2]);
patch(datenum(xp1),yp1,[.5 .5 .5],'facealpha',0.1);
title(['\fontsize{10} Most signifcant changepoint (RMS) for ',num2str(factor1),'-yr Interpolation of ' seriesname1{1}]);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
%location of most signifcant changepoint - copy paste x10^4
xp1

ax(1) = subplot(2,1,2);
plot(x3,y3ln,'ko','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',1,'LineWidth',0.2), hold on
plot(x4,y4ln,'r','LineWidth',1);
i = findchangepts(diff(y4),'Statistic','rms');
ax = gca;
xp2 = [x4(i) ax.XLim([2 2]) x4(i)];
yp2 = ax.YLim([1 1 2 2]);
patch(datenum(xp2),yp2,[.5 .5 .5],'facealpha',0.1);
title(['\fontsize{10} Most signifcant changepoint (RMS) for ',num2str(factor1),'-yr Interpolation of ' seriesname2{1}]);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
%location of most signifcant changepoint - copy paste x10^4
xp2

% Write output summary stats to file & save figure
diary('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/S1&2_Sig_CP.txt')
xp1, xp2
diary off

h51 = figure(51)
set(h51, 'Units', 'centimeters','PaperSize', [20,20]);
print(h51, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig51', '-dpdf')

%% PART 6 - WAVELET ANALYSIS

% Copyright notice for wavelet analysis
%   Copyright (C) 2002-2004, Aslak Grinsted
%
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.

%% Load the data to run independently of previous sections 
% Load up the two equally spaced time series into the matrices d1 and d2.

%seriesname={'SAm Temp An' 'Fan-Ti'};
%d1=load('data\SAm.txt');
%d2=load('data\FanTi.txt');

% Goto to section where d1 and d2 are inputs to start 

%% Or use evenly-spaced time series from before

% Using standardised detrended and PCHIP interpolated data for series 1
xq1 = x2';
vq1 = y2sln;

% Using standardised detrended and PCHIP interpolated data for series 2
xq2 = x4';
vq2 = y4sln;

%% Plotting PCHIP interpolated Z-scores time series on same figure

figure(61)
%Series 1
ax(1) = subplot(2,1,1);
plot(x1,y1ln,'o','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',1,'LineWidth',0.2), hold on
plot(xq1,vq1,'b','LineWidth',1), %hold off
xlabel(xlbl1{1});
ylabel(ylbl1{1});
yticks ([-10 -5 -2 -1 0 1 2 5 10]);
title(['\fontsize{12}',num2str(meanyr1x),'-yr IDZ (LBF) ' seriesname1{1}]);

%Series 2
ax(2) = subplot(2,1,2);
plot(x3,y3ln,'ko','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',1,'LineWidth',0.2), hold on
plot(xq2,vq2,'r','LineWidth',1), hold off
xlabel(xlbl2{1});
ylabel(ylbl2{1});
yticks ([-10 -5 -2 -1 0 1 2 5 10]);
title(['\fontsize{12}',num2str(meanyr1x),'-yr IDZ (LBF) ' seriesname2{1}]);

%link parts a to d and zoom in to show the changes
linkaxes(ax(1:2),'x')
xlim auto %([min(x1) max(x1)])
%ylim([0 max(v1-1)])
ylim auto

h61 = figure(61)
set(h61, 'Units', 'centimeters','PaperSize', [20,20]);
print(h61, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig61', '-dpdf')

%% Setting up wavelet plots
seriesname={[seriesname1{1}],[seriesname2{1}]}; %'Ln(Mn/Ti)' 'Ln(Fe/Mn)'};
d1 = [xq1,vq1];
d2 = [xq2,vq2];

%% Change timeseries into percentiles if highly bimodal - not usually
% The time series of Baltic Sea ice extent is highly bi-modal and we
% therefore transform the timeseries into a series of percentiles. The
% transformed series probably reacts 'more linearly' to climate.
% d1(:,2)=boxpdf(d1(:,2));
% d2(:,2)=boxpdf(d2(:,2));

%% Continuous wavelet transform (CWT)
% The CWT expands the time series into time
% frequency space - use standardised data input
%?The thick black contour designates the 5% significance level against 
% red noise and the cone of influence (COI) where edge effects might 
% distort the picture is shown as a lighter shade?

figure (62) %('color',[1 1 1])
tlim=[min(d1(1,1),d2(1,1)) max(d1(end,1),d2(end,1))];

subplot(2,1,1);
wt(d1);
title(['\fontsize{10}',num2str(meanyr1x),'-yr IDZ Wavelet Power Spectrum ', seriesname{1}]);
set(gca,'xlim',tlim);

subplot(2,1,2)
wt(d2)
title(['\fontsize{10}',num2str(meanyr2x),'-yr IDZ Wavelet Power Spectrum ', seriesname{2}])
set(gca,'xlim',tlim)

h62 = figure(62)
set(h62, 'Units', 'centimeters','PaperSize', [20,20]);
print(h62, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig62', '-dpdf')

%% Cross wavelet transform (XWT)
% The XWT finds regions in time frequency space where
% the time series show high common power.
% 'The 5% significance level against red noise is shown as a thick contour. 
% The relative phase relationship is shown as arrows 
% with in-phase pointing right, anti-phase pointing left' 
% where one series leads another by 90?, arrows point straight down).

figure (63) %('color',[1 1 1])
xwt(d1,d2)
title(['\fontsize{10}',num2str(factor1),'-yr IDZ Cross Wavelet Transform (XWT) ' seriesname{1} ' vs ' seriesname{2} ] )

h63 = figure(63)
set(h63, 'Units', 'centimeters','PaperSize', [20,20]);
print(h63, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig63', '-dpdf')

%% Squared Wavelet coherence (WTC)
% The WTC finds regions in time frequency space where the two
% time series co-vary (but does not necessarily have high power).

% Note: 'The area of a time frequency plot above the 5% significance level 
% is not a reliable indication of causality. Even if the scales were appropriately weighed
% for the averaging, it is possible for two series to be perfectly
% correlated at one specific scale while the area of significant
% correlation is much less than 5%. However, if the significant region 
% in the plot is very extensive at/across a specific timescale, 
% it is very unlikely that this is simply by chance' 
% from: https://www.glaciology.net/pdf/Grinsted-npg2004-wavelet-coherence.pdf 

% Advantages (also from above)
% 'the WTC can be thought of as the local correlation 
% between the time series in time frequency space. 
% XWT reveals high common power, WTC finds locally phase locked behavior. 
% The more desirable features of the WTC come at the price of being slightly
% less localized in time frequency space. 
% The significance level of the WTC has to be determined using Monte Carlo
% method'
% from: https://www.glaciology.net/pdf/Grinsted-npg2004-wavelet-coherence.pdf

figure(64) %('color',[1 1 1])
%ax(1) = subplot(3,1,1);
wtc(d1,d2)
title(['\fontsize{10}',num2str(factor1),'-yr IDZ Wavelet Transform Coherence (WTC) ' seriesname{1} ' vs ' seriesname{2} ] );

h64 = figure(64)
set(h64, 'Units', 'centimeters','PaperSize', [20,20]);
print(h64, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig64', '-dpdf')

%% Wavelet coherence (WTC) plotted agianst original/interpolated time series
% The WTC finds regions in time frequency space where the two
% time series co-vary (but does not necessarily have high power).

figure(65) %('color',[1 1 1])
ax(1) = subplot(3,1,1);
wtc(d1,d2)
title(['\fontsize{10}',num2str(factor1),'-yr IDZ Wavelet Transform Coherence (WTC) ' seriesname{1} ' vs ' seriesname{2} ] );

ax(2) = subplot(3,1,2);
plot(x1,y1ln,'o','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',1,'LineWidth',0.2), hold on
plot(xq1,vq1,'b','LineWidth',1), %hold off
ylabel(ylbl1{1});

ax(3) = subplot(3,1,3);
plot(x3,y3ln,'o','MarkerEdgeColor','#a9a9a9','MarkerFaceColor','#a9a9a9','MarkerSize',1,'LineWidth',0.2), hold on
plot(xq2,vq2,'r','LineWidth',1), hold off
xlabel(xlbl2{1});
ylabel(ylbl2{1});

%link parts a to d and zoom in to show the changes
linkaxes(ax(1:3),'x')
xlim auto
ylim auto

h65 = figure(65)
set(h65, 'Units', 'centimeters','PaperSize', [20,20]);
print(h65, '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/Fig65', '-dpdf')

%% Write all figures to folder

FolderName = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp'  % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  savefig(fullfile(FolderName, [FigName '.fig']));
end

%% Compile output files in order of final output file
output1 = [series1x, series1ynorm, series1ydl, series1ydlnorm];
output2 = [series2x, series2ynorm, series2ydl, series2ydlnorm];
output3 = [x1reshape, y1reshape, y1normreshape, x2reshape, y2normreshape, y2s, y2sl, y2sln, y4s, y4sl, y4sln];

out1 = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/output1.csv';
xlswrite(out1, output1);
out2 = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/output2.csv';
xlswrite(out2, output2);
out3 = '/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Temp/output3.csv';
xlswrite(out3, output3);

%% close all figure windows and end

close all
 