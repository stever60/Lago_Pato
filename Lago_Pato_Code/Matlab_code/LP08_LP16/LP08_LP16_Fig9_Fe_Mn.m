%% Clear everything and start again
clc
clear all
close all
 
%%  Paths for Macbook
addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/Wavelet/grinsted-wavelet-coherence-d987ea4');
addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/Wavelet/grinsted-wavelet-coherence-d987ea4/data');
addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis');
addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Data')
addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/fLOESS')

%% Paths for Windows desktop
% Add paths for functions and data needed in this sheet on Windows 10 desktop computer
% addpath('F:\Dropbox\BAS\Data\Matlab\Matlab_projects\Wavelet\grinsted-wavelet-coherence-d987ea4');
% addpath('F:\Dropbox\BAS\Data\Matlab\Matlab_projects\Wavelet\grinsted-wavelet-coherence-d987ea4\data');
% addpath('F:\Dropbox\BAS\Data\Matlab\Matlab_projects\PeakAnalysis');
% addpath('F:\Dropbox\BAS\Data\Matlab\Matlab_projects\PeakAnalysis\Data')
% addpath('F:\Dropbox\BAS\Data\MATLAB\Matlab_projects\fLOESS')

%% Import and check 2 data series
%Need to insert import sequences from dataplots.m and correct seriesname
%and y axis label (ylbl) for plots
%load datafiles into matlab
%remove NaNs code
%out = A(all(~isnan(A),2),:); % for nan - rows
%out = A(:,all(~isnan(A)));   % for nan - columns

%Series 1 and remove NaN rows
load LP08_Fig9_LnFe_Mn.txt
series1 = LP08_Fig9_LnFe_Mn(all(~isnan(LP08_Fig9_LnFe_Mn),2),:);
%Series 2
load LP16_Fig9_LnFe_Mn.txt
series2 = LP16_Fig9_Fe_Mn(all(~isnan(LP16_Fig9_Fe_Mn),2),:);

%% Setting up two series
% Create common file name and take natural log of the data
% y = log(x) returns the natural logarithm ln(x) of each element in the array

%Series 1
format long;
series1x = (series1(:,1));
series1y = (series1(:,2));
%series1y = log(LP08_I_Mn_Ti(:,2));
x1 = series1x;
y01 = series1y; %original data series 1

%check output - if needed - ; stops the output displaying in command window
series1x
series1y

% Series 2
format long;
series2x = (series2(:,1));
series2y = (series2(:,2));
%series2y = log(LP08_I_Fe_Mn(:,2));
x3 = series2x;
y02 = series2y; %original data series 2

% Create series labels and x, y axis labels for all plots
seriesname1={'LP08 Ln(Fe/Mn)'};
ylbl1 = {'LP08: Ln(Fe/Mn)'};
xlbl1 = {'Age (cal yr BP)'};

seriesname2={'LP16 Ln(Fe/Mn)'};
ylbl2 = {'LP16: Ln(Fe/Mn)'};
xlbl2 = {'Age (cal yr BP)'};

% special label if series 1 and seris 2 are use in same figure
ylblspecial1 = {'LP08: Ln(Fe/Mn) (blue), LP16: Ln(Fe/Mn) (orange)'};

%% Standardize and centre y (Z-scores) 
% subtracting from mean and dividing by stdev for series 1 and 2
series1ynorm = (series1y-mean(series1y)) %/std(series1y)
series2ynorm = (series2y-mean(series2y)) %/std(series2y)

%write to file 
s1norm = 's1n.csv';
xlswrite(s1norm,series1ynorm);
s1xnorm = 's1xn.csv';
xlswrite(s1xnorm,series1x);
s2norm = 's2n.csv';
xlswrite(s2norm,series2ynorm);
s2xnorm = 's2xn.csv';
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
legend('Measured', 'Standardized');
legend('Location','northeast')
legend ('boxon');

ax(2) = subplot(2,1,2);
plot(series2x,series2y,series2x,series2ynorm);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
title(['\fontsize{12}','Time series plot for ' seriesname2{1}]);
xlim([min(series2x) max(series2x)]);
xlim auto;
ylim auto;
legend('Measured', 'Standardized');
legend('Location','northwest')
legend ('boxon');

%link parts and zoom in to show the changes
linkaxes(ax(1:2),'x');
%xlim([min(series1x) max(series1x)])
xlim ([min(-100) max(20000)])
xlim auto
ylim auto

%% Establish min and max limits for x-axis (Years) and interval spacing 

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

% smooth to, for example, nearest 10 years from one year rounded data 
% need to enter factor = 10 as data input here
factor1 = meanyr1*(10/meanyr1)
factor2 = meanyr2*(10/meanyr2)
meanyr1x = factor1
meanyr2x = factor2

%1yr dataset to make two vectors of the same size
factor3 = meanyr1*(1/meanyr1)
factor4 = meanyr2*(1/meanyr2)
meanyr3x = factor3
meanyr4x = factor4

%% Plot time interval summary and add values
figure(2)
ax(1) = subplot(2,1,1);
plot(intv1)
hold on
plot(intv1sm)
title(['\fontsize{12}','Time interval plot for ' seriesname1{1}]);
xlabel('Sample no.');
ylabel(ylbl1{1});
ylim ([min(0) max(20)]);
legend(['Mean interval (years±2s) = ',num2str(meanyr1),  '±',num2str(stdevyr1)],['Mov. Ave. mean = ',num2str(meanyr1sm),  ' years']);
legend('Location','northeast')
legend ('boxon');
hold off

ax(2) = subplot(2,1,2);
plot(intv2)
hold on
plot(intv2sm)
title(['\fontsize{12}','Time interval plot for ' seriesname2{1}]);
xlabel('Sample no.');
ylabel(ylbl2{1});
ylim ([min(0) max(20)]);
legend(['Mean interval (years±2s) = ',num2str(meanyr2),  '±',num2str(stdevyr2)],['Mov. Ave. mean = ',num2str(meanyr2sm),  ' years']);
legend('Location','northeast')
legend ('boxon');
hold off

%link parts and zoom in to show the changes
linkaxes(ax(1:2),'x');
%xlim([min(series1x) max(series1x)])
xlim ([min(-100) max(20000)])
xlim auto
ylim auto

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

x1reshape = reshape(t1,[984,1]);
x2reshape = reshape(t2,[1303,1]);

y1reshape = reshape(s1inter1,[984,1]);
y2reshape = reshape(s2inter1,[1303,1]);
y1normreshape = reshape(s1norminter,[984,1]);
y2normreshape = reshape(s2norminter,[1303,1]);

%write to file 
sx1inter = 'sx1inter.csv';
xlswrite(sx1inter,x1reshape);
sx2inter = 'sx2inter.csv';
xlswrite(sx2inter,x2reshape);

sy1inter = 'sy1inter.csv';
xlswrite(sy1inter,y1reshape);
sy2inter = 'sy2inter.csv';
xlswrite(sy2inter,y2reshape);

sy1norminter = 'sy1norminter.csv';
xlswrite(sy1norminter,y1normreshape);
sy2inter = 'sy2norminter.csv';
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


%% Plot the normalised and interpolated data 
% plot series 1 data with linear (red line plot) and then pchip (blue line plot) interpolations
% with datapoints plotted as black (k) filled circles (o)

figure
%plot series 2 data with linear and pchip interpolations
ax(1) = subplot(2,1,1);
plot(x3,y02n,'k.','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',2,'LineWidth',0.2), hold on
plot(x4,y2ni, 'r-','LineWidth',1), hold off
xlabel(xlbl2{1});
ylabel(ylbl2{1});
title(['\fontsize{12}',num2str(meanyr2x),'-yr PCHIP Interpolated Z-scores for  ' seriesname2{1}]);
xlim auto
ylim ([min(-6) max(4)])
%ylim auto %([min(-1) max(1)]);

ax(2) = subplot(2,1,2);
plot(x1,y01n,'.','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',2,'LineWidth',0.2), hold on
plot(x2,y1ni, 'b-','LineWidth',1), hold off
xlabel(xlbl1{1});
ylabel(ylbl1{1});
title(['\fontsize{12}',num2str(meanyr1x),'-yr PCHIP Interpolated Z-scores for ' seriesname1{1}]);
xlim auto
ylim ([min(-6) max(4)])
%ylim auto %([min(-1) max(1)]);

%link parts and zoom in to show the changes
linkaxes(ax(1:2),'x');
%xlim([min(series1x) max(series1x)])
xlim ([min(-100) max(20000)])
ylim ([min(-6) max(4)])
%xlim auto
%ylim auto

%% Detrend the y data 

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

detrended1y = 'series1ydl.csv';
xlswrite(detrended1y,series1ydl);

detrended2y = 'series2ydl.csv';
xlswrite(detrended2y,series2ydl);

detrended1ynorm = 'series1ydlnorm.csv';
xlswrite(detrended1ynorm,series1ydlnorm);

detrended2ynorm = 'series2ydlnorm.csv';
xlswrite(detrended2ynorm,series2ydlnorm);

%% Plot detrended data to check 
% series1yd %this is long but can copy from data window if needed elsewhere
figure(3)
%Linking axis 
ax(1) = subplot(4,1,1);
plot(series1x,series1y,series1x,series1yd,series1x,series1ydl)
hold on
%plot(x1,y1,'.'); %check these are the same by running this as well
title(['\fontsize{12}','Detrended time series plot for ' seriesname1{1}]);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
legend('Measured','Detrended (mean)','Detrended (LBT)');
legend('Location','northwest')
legend ('boxon');

ax(2) = subplot(4,1,2);
plot(series2x,series2y,series2x,series2yd,series2x,series2ydl)
hold on
%plot(x3,y3,'.'); %check these are the same by running this as well
title(['\fontsize{12}','Detrended time series plot for ' seriesname2{1}]);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
legend('Measured','Detrended (mean)','Detrended (LBT)');
legend('Location','northwest')
legend ('boxon');

ax(3) = subplot(4,1,3);
plot(series1x,series1ynorm,series1x,series1ydnorm,series1x,series1ydlnorm)
hold on
%plot(x1,y1,'.'); %check these are the same by running this as well
title(['\fontsize{12}','Detrended time series plot of Z-scores for ' seriesname1{1}]);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
legend('Z-scores','Z-scores Detrend (mean)','Z-scores Detrend (LBT)');
legend('Location','northwest')
legend ('boxon');

ax(4) = subplot(4,1,4);
plot(series2x,series2ynorm,series2x,series2ydnorm,series2x,series2ydlnorm)
hold on
%plot(x1,y1,'.'); %check these are the same by running this as well
title(['\fontsize{12}','Detrended time series plot of Z-scores for ' seriesname2{1}]);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
legend('Z-scores','Z-scores Detrend (mean)','Z-scores Detrend (LBT)');
legend('Location','northwest')
legend ('boxon');

%link parts and zoom in to show the changes
linkaxes(ax(1:4),'x');
%xlim([min(series1x) max(series1x)])
xlim auto
ylim auto

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
y2 = series1p; %detrended and factor interpolated data series 1
y2l = series1pl; %detrended data series 1 - subtracted from linear best fit
y4 = series2p; %detrended and factor interpolated data series 2 - subtracted from mean
y4l = series2pl; %detrended and factor interpolated data series 2 - subtracted from linear

%Standardized, Detrended & interpolated data
y2ln = series1plnorm; %std and detrended data series 1 - subtracted from linear best fit
y4ln = series2plnorm; %stdized and detrended and factor interpolated data series 2 - subtracted from linear

%reshaped interpolated detrended y series into a vector with n rows x 1
%column from n columns x 1 row - to use later on for regression etc.

%Series 1
x2s = reshape(x2,[984,1]);
y2s = reshape(y2,[984,1]);
y2sl = reshape(y2l,[984,1]);
y2sln = reshape(y2ln,[984,1]);

%Series 2
x4s = reshape(x4,[1303,1]);
y4s = reshape(y4,[1303,1]);
y4sl = reshape(y4l,[1303,1]);
y4sln = reshape(y4ln,[1303,1]);

%Nomalise outputs - if needed
%y2slnorm = (y2sl-mean(y2sl))/std(y2sl);
%y4slnorm = (y4sl-mean(y4sl))/std(y4sl);

%write these values to csv - example - or copy paste from command window 
%need to open a blank excel spreadsheet first for this to work

%series 1
interx2s = 'x2s.csv';
xlswrite(interx2s,x2s);

intery2s = 'y2s.csv';
xlswrite(intery2s,y2s);
intery2sl = 'y2sl.csv';
xlswrite(intery2sl,y2sl);
intery2sln = 'y2sln.csv';
xlswrite(intery2sln,y2sln);
%interpolatedy2slnorm = 'y2slnorm.csv';
%xlswrite(interpolatedy2slnorm,y2slnorm);

%series 2
interx4s = 'x4s.csv';
xlswrite(interx4s,x4s);

intery4s = 'y4s.csv';
xlswrite(intery4s,y4s);
intery4sl = 'y4sl.csv';
xlswrite(intery4sl,y4sl);
intery4sln = 'y4sln.csv';
xlswrite(intery4sln,y4sln);
%interpolatedy4slnorm = 'y4slnorm.csv';
%xlswrite(interpolatedy4slnorm,y4slnorm);

%% Plot the detrended and interpolated data 

% plot series 1 data with linear (red line plot) and then pchip (blue line plot) interpolations
% with datapoints plotted as black (k) filled circles (o)

figure (4)
ax(1) = subplot(4,1,1);
plot(x1,y1l,'.','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',2,'LineWidth',0.2), hold on
plot(x2,y2l, 'r-','LineWidth',1), hold off
%plot(x2,y2lst, 'r-','LineWidth',1), hold on
xlabel(xlbl1{1});
ylabel(ylbl1{1});
title(['\fontsize{12}',num2str(meanyr1x),'-yr PCHIP Interpolated/Detrended (LBT) for ' seriesname1{1}]);
xlim auto
ylim auto %([min(-1) max(1)]);

%plot series 2 data with linear and pchip interpolations
ax(2) = subplot(4,1,2);
plot(x3,y3l,'k.','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',2,'LineWidth',0.2), hold on
plot(x4,y4l, 'r-','LineWidth',1), hold off
xlabel(xlbl2{1});
ylabel(ylbl2{1});
title(['\fontsize{12}',num2str(meanyr2x),'-yr PCHIP Interpolated/Detrended for ' seriesname2{1}]);
xlim auto
ylim auto %([min(-1) max(1)]);

%plot series 1 and 2 standardized data with linear pchip interpolation
ax(3) = subplot(4,1,3);
plot(x1,y1ln,'.','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',2,'LineWidth',0.2), hold on
plot(x2,y2ln, 'b-','LineWidth',1), hold off
xlabel(xlbl1{1});
ylabel(ylbl1{1});
title(['\fontsize{12}',num2str(meanyr1x),'-yr PCHIP Interpolated/Detrended z-scores (LBT) for ' seriesname1{1}]);
xlim auto
ylim auto %([min(-1) max(1)]);

ax(4) = subplot(4,1,4);
plot(x3,y3ln,'k.','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',2,'LineWidth',0.2), hold on
plot(x4,y4ln, 'b-','LineWidth',1), hold off
xlabel(xlbl2{1});
ylabel(ylbl2{1});
title(['\fontsize{12}',num2str(meanyr2x),'-yr PCHIP Interpolated/Detrended z-scores (LBT) for ' seriesname2{1}]);
xlim auto
ylim auto %([min(-1) max(1)]);

%% Scatterplot of interpolated-year detrended data compared with 10 or 100 year interpolated data 
figure (5)
ax(1) = subplot(2,1,1)
%scatter(series1pl1,series2pl1, '.','MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',0.2);
hold on;
scatter(y2l,y4l, 'o','MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',0.2);xlabel(ylbl1{1});
hold off;
ylabel(ylbl2{1});
title(['\fontsize{12}', seriesname1{1}, ' vs ',seriesname2{1}]);
legend('10-year detrended (LBT) PCHIP Interpolated');
legend('Location','northeast')
legend ('boxoff');
xlim auto
ylim auto

ax(2) = subplot(2,1,2)
%scatter(series1plnorm,series2plnorm, '.','MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',0.2);
hold on;
scatter(y2ln,y4ln, 'o','MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',0.2);xlabel(ylbl1{1});
hold off;
ylabel(ylbl2{1});
title(['\fontsize{12}', seriesname1{1}, ' vs ',seriesname2{1},'PCHIP,Std & detrended']);
%legend('10-year standardised & detrended (LBT) PCHIP Interpolated');
%legend('Location','northeast')
%legend ('boxoff');
xlim auto
ylim auto

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
figure(6)
ax(1) = subplot(2,2,1);
[Pxx1,f1] = periodogram(y2ln,[],256,1/meanyr1x);
plot(f1,Pxx1,'k');
xlim ([min(0) max(0.04)]) %check whole range using 'xlim auto' first and then rescale
ylim ([min(0) max(200)]) %check whole range using 'ylim auto' first and then rescale
xlabel('Frequency')
ylabel('Power')
title(['\fontsize{10}','Auto-spectrum of ' seriesname1{1}]);

%series 1 - interpolated years plot 
ax(2) = subplot(2,2,2);
[Pxx1,f1] = periodogram(y2ln,[],256,1/meanyr1x);
semilogx(1./f1,Pxx1,'k');
xlim ([min(0) max(10000)])
xticks auto %use auto first to check range
ylim ([min(0) max(200)]) %use auto first to check range
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
plot(1./f1(locsf1),pksf1+2,'v','markerfacecolor','w','markeredgecolor','b', 'MarkerSize',4);
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
ylim ([min(0) max(200)]) %check whole range using 'ylim auto' first
xlabel('Frequency')
ylabel('Power')
title(['\fontsize{10}','Auto-spectrum of ' seriesname2{1}]);

%series 2 - interpolated years plot
ax(4) = subplot(2,2,4);
[Pxx2,f2] = periodogram(y4ln,[],248,1/meanyr2x);
semilogx(1./f2,abs(Pxx2),'k')
xlim ([min(0) max(10000)]) %use auto first to check range
xticks auto
ylim ([min(0) max(200)]) %use auto first to check range
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
plot(1./f2(locsf2),pksf2+0.2,'v','markerfacecolor','w','markeredgecolor','b', 'MarkerSize',4);
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
s1period = 's1period.csv';
xlswrite(s1period,period1);

s2period = 's2period.csv';
xlswrite(s2period,period2);

%% Axis formatting example - dont run 
% ax = gca
% ax.TickDir = 'out'
% ax.FontSize = 10;
% ax.TickLength = [0.04 0.04]

%% Magnitude squared coherence of the time series (MSC)
% use this to identify frequency domain correlations between two time series
% phase estimates in the cross-spectrum are where significant frequency-
% domain correlation exists - mcohere (hamming type) uses Welch's overlapped
% segment averaging (WOSA) and multitaper techniques wihth a window length
% of 100 samples, equivalent to 10 periods of a 100 Hz signal and overlap
% of 80 - therefore need to mulitply the meanyr interpolated interval by 10
% to give an output that is in Hz

figure (7)
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
[Cxy,f4]= mscohere(y2ln,y4ln,hamming(100),80,meanyr1x*(factor1)^2);
% x-axis in years
xf4 = (1./f4)*factor1^2;
semilogx(xf4,Cxy,'k');
xlabel('Years')
ylabel('Magnitude Squared Coherence')
title(['\fontsize{8}','Coherence of ' seriesname1{1}, ' and ',seriesname2{1}]);
xlim ([min(0) max(10000)]);
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
plot(xf4(locsf4),pksf4+0.1,'v','markerfacecolor','w','markeredgecolor','b', 'MarkerSize',4);
hold off;
%location of coherent peaks in years and peak value
mscpeak = xf4(locsf4)
pksf4
%link y axes for Series 1 plots
linkaxes([ax(1),ax(2)],'y')

%Write peak locations to file 
mscpk = 'msc_peak_years.csv';
xlswrite(mscpk,mscpeak);

% Frequency cross spectrum 
ax(3) = subplot(2,2,3);
[Pxy,f3]= cpsd(series1pl,series2pl,hamming(100),80,meanyr1x*10);
plot(f3,abs(Pxy),'k');
xlabel('Frequency (Hz)')
ylabel('Power')
xlim auto %([min(0) max(0.1)])
ylim auto %([min(0) max(1.2)])
title(['\fontsize{8}','Cross-Spectrum of ' seriesname1{1}, ' & ',seriesname2{1}]);
hold on

ax(4) = subplot(2,2,4);
[Pxy,f3]= cpsd(series1pl,series2pl,hamming(100),80,meanyr1x*10);
%create new axis from factor^2 = 100
xf3 = (1./f3)*factor1^2; 
semilogx(xf3,abs(Pxy),'k');
xlabel('Years')
ylabel('Power')
title(['\fontsize{8}','Cross-Spectrum of ' seriesname1{1}, ' & ',seriesname2{1}]);
xlim ([min(0) max(10000)]);
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
plot(xf3(locsf3),pksf3+0.00005,'v','markerfacecolor','w','markeredgecolor','b', 'MarkerSize',4);
hold off;
%location of cross spectrum peaks in years and peak value
xf3(locsf3)
pksf3

%Write peak locations to file 
cspect = 'cspect.csv';
xlswrite(cspect,xf3(locsf3));

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

%% Phase lag (in radians) of cross spectrum
% plot is of the phase of the cross spectrum
% which indicates frequencies with sig coherence between two series
% and shows known phase between their sinusoidal components of input
% signals

figure(8)
phase = angle(Pxy);
plot(f3,-angle(Pxy)/pi);
xlabel('Frequency (Hz)')
ylabel('Lag \times\pi rad')
title(['\fontsize{12}','Cross spectrum Phase of ' seriesname1{1}, ' & ',seriesname2{1}]);
xlim auto %([min(0) max(0.1)])
ylim auto %([min(-4) max(4)])

%% Evolutionary power spectrum - spectrogram map
%Computes a short-time Forier transform (STFT) of overlapping segments of the
%time series using a Hamming window of 400 datapoints, 50 data point overlap 
% STFT 256 resolution (higher produces a smoother plot) and frequency set
% to sample as above 
% not sure if this is correct - need better interpretation of what this
% means - check axes are labelled properly 

figure(9)
spectrogram(series1pl,64,50,256,1/factor1);
title(['\fontsize{12}','Evolutionary Power Spectrum of ' seriesname1{1}]);
xlabel('Frequency (1/yr)')
ylabel('Time (yr)')
colormap jet;
xlim auto %([min(0) max(150)]);
ylim auto %([min max]);

figure(10)
spectrogram(series2pl,64,50,256,1/factor1);
title(['\fontsize{12}','Evolutionary Power Spectrum of ' seriesname2{1}]);
xlabel('Frequency (1/yr)')
ylabel('Time (yr)')
colormap jet;
xlim auto %([min(0) max(150)]);
ylim auto %([min max]);

%% Least squares polynomial fitting
%Plot series 1 data and compare to interpolated plot

figure(20)
ax(1) = subplot(4,1,1);
plot(x1,series1ydl,'.','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',2,'LineWidth',0.2), hold on
plot(x2,y2l,'r','LineWidth',1), hold off
xlabel(xlbl1{1});
ylabel(ylbl1{1});
title(['\fontsize{12}',num2str(meanyr1x),'-yr PCHIP Interpolated/Detrended (LBT) of' seriesname1{1}]);
axis auto %([yrmin1 yrmax1 min(y1) max(y1)]);
ylim auto %([min(-1) max(1)]);
xlim auto %([min(x1) max(x1)]);
%yticks ([-10 -5 -2 -1 0 1 2 5 10])
 
ax(2) = subplot(4,1,2);
plot(x3,series2ydl,'.','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',2,'LineWidth',0.2), hold on
plot(x4,y4l,'b','LineWidth',1), hold off
xlabel(xlbl2{1});
ylabel(ylbl2{1});
title(['\fontsize{12}',num2str(meanyr2x),'-yr PCHIP Interpolated/Detrended (LBT) of' seriesname2{1}]);
axis([yrmin1 yrmax1 min(y1) max(y1)]);
xlim auto
ylim auto %([min(-1) max(1)]);
xlim auto %([min(x1) max(x1)]);
%yticks ([-10 -5 -2 -1 0 1 2 5 10])

ax(3) = subplot(4,1,3);
plot(x1,series1ydlnorm,'.','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',2,'LineWidth',0.2), hold on
plot(x2,y2ln,'r','LineWidth',1), hold off
xlabel(xlbl1{1});
ylabel(ylbl1{1});
title(['\fontsize{12}',num2str(meanyr1x),'-yr PCHIP Interpolated/Detrended (LBT) z-scores of' seriesname1{1}]);
axis auto %([yrmin1 yrmax1 min(y1) max(y1)]);
ylim auto %([min(-1) max(1)]);
xlim auto %([min(x1) max(x1)]);
yticks ([-10 -5 -2 -1 0 1 2 5 10]) %auto %([0:100:1000]);
 
ax(4) = subplot(4,1,4);
plot(x3,series2ydlnorm,'.','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',2,'LineWidth',0.2), hold on
plot(x4,y4ln,'b','LineWidth',1), hold off
xlabel(xlbl2{1});
ylabel(ylbl2{1});
title(['\fontsize{12}',num2str(meanyr2x),'-yr PCHIP Interpolated/Detrended (LBT) z-scores of' seriesname2{1}]);
axis([yrmin1 yrmax1 min(y1) max(y1)]);
xlim auto
ylim auto %([min(-1) max(1)]);
xlim auto %([min(x1) max(x1)]);
yticks ([-10 -5 -2 -1 0 1 2 5 10])

linkaxes([ax(1),ax(2), ax(3), ax(4)], 'x');
linkaxes([ax(1),ax(2), ax(3), ax(4)], 'y');

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

figure(21)
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

title(['\fontsize{10}Least squares polynomial fit for ', seriesname1{1}]);
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
title(['\fontsize{10}Least squares polynomial fit for ',num2str(meanyr1x),'-yr PCHIP Interpolation of Detrended' seriesname1{1}]); 
plot(x2,vq2fit,'linewidth',1), hold off;
xlim auto
ylim auto

%reshape the interpolated detrended y series into a vector with n rows x 1
%column from n columns x 1 row - to use later on for regression etc.
x2s = reshape(x2,[984,1])
y2s = reshape(y2,[984,1])
x4s = reshape(x4,[1303,1])
y4s = reshape(y4,[1303,1])

x2s = reshape(x2,[984,1])
y2sl = reshape(y2l,[984,1])
x4s = reshape(x4,[1303,1])
y4sl = reshape(y4l,[1303,1])

% regress and plot the data
ax(4) = subplot(4,1,4);
residuals2 = y2s - vq2fit;
[b2,bint2,r2] = regress(y2s,x2s);
plot(x2s,residuals2,'linewidth',1)%residuals from fit above
xlabel(xlbl1{1});
ylabel('Residuals');
xlim auto
ylim auto
format long
[p2,DW2] = dwtest(r2,x2s) %exact p values in two tailed test
title(['\fontsize{10}Residuals p=',num2str(p2),'; DW=',num2str(DW2)])

%% Summary stats output for Series 1
% write these values to file
%filenamestats1 = 'Polyfit_stats_Series1.csv';
%series1stats = {'p-value','DW stat'; p2, DW2}; 
%csvwrite(filenamestats1, series1stats);

%% Polyfit Series 2
%Series 2 plots 

figure(22)
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
[b4,bint4,r4] = regress(y4s,x4s);
plot(x4,residuals4,'linewidth',1)%residuals from fit above
xlabel(xlbl2{1});
ylabel('Residuals');
xlim auto
ylim auto
format long
[p4,DW4] = dwtest(r4,x4s) %exact p values in two tailed test
title(['\fontsize{10}Residuals p=',num2str(p4),'; DW=',num2str(DW4)])

linkaxes([ax(1),ax(2),ax(3),ax(4)],'x');

%% Summary stats output for Series 2
% write these values to file
% filenamestats4 = 'Polyfit_stats_Series1.csv';
% series2stats = {'p-value','DW stat'; p4, DW4}; 
% csvwrite(filenamestats4, series2stats);

%% Serial Autocorrelation
% autocorrelation using residual plots
% using 'residuals' from polynomial trend fitting plot

figure(25)
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
title(['\fontsize{10}AC of Std-Detrend-Int ', seriesname1{1}]);

%plot the autocorrelation for derended and interpolated series 2
ax(4) = subplot(2,2,4);
acf4 = autocorr(y4ln);
%x = autocorr(v1) for original data
semilogx(acf4);
xlabel('Lag time (years)');
ylabel('Correlation coefficient');
title(['\fontsize{10}AC of Std-Detrend-Int ', seriesname2{1}]);


%% Peak identification

% Series 1
figure(27)
% This is a 4-part figure

% figure a - peak identification to find all the data points as the peaks
ax(1) = subplot(4,1,1);
%plot(x1,y1); hold on;
plot(x2,y2sln); hold on;
[pks0,locs0] = findpeaks(y2sln,'MinPeakProminence',0);
% offset values of peak heights for plotting
plot(x2(locs0),pks0+0.4,'bv','markerfacecolor',[0 0 1], 'MarkerSize',3);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
locs0yrs= locs0*factor1
locs0yrs
locs0mean = round(mean(diff(locs0yrs)))
title(['\fontsize{10} All peaks (',num2str(meanyr1x),'-yr IDZ data): ',seriesname1{1},' Mean peak diff = ',num2str(locs0mean),' yrs'])


% figure b
% 1s filtering and caluclating distance between peaks in first plot
% sets MinPeakProminence value to >stdev of interpolated data interval
MPP = std(intv1) 
ax(2) = subplot(4,1,2);
%plot(x1,y1); hold on;
plot(x2,y2sln); hold on;
[pks1,locs1] = findpeaks(y2sln,'MinPeakProminence',MPP);
% calculate offset values of peak heights for plotting
plot(x2(locs1),pks1+0.4,'bv','markerfacecolor',[0 1 1], 'MarkerSize',3);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
locs1yrs= locs1*factor1
locs1yrs
locs1mean = round(mean(diff(locs1yrs)))
title(['\fontsize{10} Prominent IDZ Peaks [Trough-Peak Diff >',num2str(2*MPP),' (1s)] Mean peak diff = ',num2str(locs1mean),' yrs'])%adds variable to title - remember []


% figure c 
% 2s filtering and caluclating distance between peaks in first plot
ax(3) = subplot(4,1,3);
[pks2,locs2] = findpeaks(y2sln,'MinPeakProminence',2*MPP);
%plot(x1,y1); hold on;
plot(x2,y2sln); hold on;
% offset values of peak heights for plotting
plot(x2(locs2),pks2+0.4,'bv','markerfacecolor',[1 1 0], 'MarkerSize',3);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
locs2yrs= locs2*factor1
locs2yrs
locs2mean = round(mean(diff(locs2yrs)))
title(['\fontsize{10} IDZ peak-trough variance >2s Mean peak diff = ',num2str(locs2mean),' yrs']) %adds variable to title

% figure d 
% detecting saturated peaks using 2-sigma as a vertical constraint 
% versus T1 interpolation as a horizontal peak-width constraint
ax(4) = subplot(4,1,4);
%plot(x1,y1); hold on;
plot(x2,y2sln); hold on;
TH = (MPP)/T1 %sets threshold value to stdev / interpolation years e.g., 10 or 50 years
[pks3,locs3] = findpeaks(y2sln,'threshold',TH);
% offset values of peak heights for plotting
plot(x2(locs3),pks3+0.4,'bv','markerfacecolor',[1 1 0], 'MarkerSize',3);
peakInterval3 = diff(locs3);
xlabel(xlbl1{1});
ylabel(ylbl1{1});
locs3yrs= locs3*factor1
locs3yrs
locs3mean = round(mean(diff(locs3yrs)))
title(['\fontsize{10} IDZ w/o saturated peaks: Threshold > ',num2str(TH), ' (1s/int.) Mean peak diff = ',num2str(locs3mean),' yrs'])%adds variable to title


%link parts a to d and zoom in to show the changes
linkaxes(ax(1:4),'xy')
xlim([min(x1) max(x1)])
%ylim([0 max(v1-1)])
ylim auto

%% Series 2
figure(28)
% This is a 4-part figure 
% figure a - peak identification to find all the data points as the peaks
ax(1) = subplot(4,1,1);
%plot(x3,y3); hold on;
plot(x4,y4sln); hold on;
[pks4,locs4] = findpeaks(y4sln,'MinPeakProminence',0);
% offset values of peak heights for plotting
plot(x4(locs4),pks4+0.5,'bv','markerfacecolor',[0 0 1], 'MarkerSize',3);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
locs4yrs= locs4*factor1
locs4yrs
locs4mean = round(mean(diff(locs4yrs)))
title(['\fontsize{10} All peaks (',num2str(meanyr2x),'-yr IDZ data)',seriesname2{1},' Mean peak diff = ',num2str(locs4mean),' yrs'])

% figure b
% 1s filtering and caluclating distance between peaks in first plot
% sets MinPeakProminence value to >stdev of interpolated data
MPP = std(intv2) %stdev_all as alternative
ax(2) = subplot(4,1,2);
%plot(x3,y3); hold on;
plot(x4,y4sln); hold on;
[pks5,locs5] = findpeaks(y4sln,'MinPeakProminence',MPP);
% calculate offset values of peak heights for plotting
plot(x4(locs5),pks5+0.5,'bv','markerfacecolor',[0 1 1], 'MarkerSize',3);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
locs5yrs= locs5*factor1
locs5yrs
locs5mean = round(mean(diff(locs5yrs)))
title(['\fontsize{10} Prominent IDZ Peaks [Trough-Peak Difference >',num2str(2*MPP),' (1s)] Mean peak diff = ',num2str(locs5mean),' yrs'])%adds variable to title - remember []
 
% figure c 
% 2s filtering and caluclating distance between peaks in first plot
ax(3) = subplot(4,1,3);
[pks6,locs6] = findpeaks(y4sln,'MinPeakProminence',2*MPP);
%plot(x3,y3); hold on;
plot(x4,y4sln); hold on;
% offset values of peak heights for plotting
plot(x4(locs6),pks6+0.5,'bv','markerfacecolor',[1 1 0], 'MarkerSize',3);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
locs6yrs= locs6*factor1
locs6yrs
locs6mean = round(mean(diff(locs6yrs)))
title(['\fontsize{10} IDZ peak-trough variance >2s Mean peak diff = ',num2str(locs6mean),' yrs'])%adds variable to title

% figure d 
% detecting saturated peaks using 2-sigma as a vertical constraint 
% versus T1 interpolation as a horizontal peak-width constraint
ax(4) = subplot(4,1,4);
%plot(x3,y3); hold on;
plot(x4,y4sln); hold on;
TH = (MPP)/T2 %sets threshold value to 2xstdev / interpolation value e.g., 10 or 50 years
[pks7,locs7] = findpeaks(y4sln,'threshold',TH);
% offset values of peak heights for plotting
plot(x4(locs7),pks7+0.5,'bv','markerfacecolor',[1 1 0], 'MarkerSize',3);
peakInterval4 = diff(locs7);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
locs7yrs= locs7*factor1
locs7yrs
locs7mean = round(mean(diff(locs7yrs)))
title(['\fontsize{10} IDZ w/o saturated peaks: Threshold > ',num2str(TH), ' (1s/int.) Mean peak diff = ',num2str(locs7mean),' yrs'])%adds variable to title
 
%link parts a to d and zoom in to show the changes
linkaxes(ax(1:4),'xy')
xlim([min(x2) max(x2)])
%ylim([0 max(v1-1)])
ylim auto
 
%% PART 4 PEAK FILTERING

%Peak-prominence >mean+2s of detrended and interpolated data 
figure (41)
ax(1) = subplot(2,1,1);
%series 1
%measuring width of the peaks using half prominence as a reference
P1 = meanyr1+(2*stdevyr1); %sets MinPeakProminence value to >stdev
findpeaks(y2,x2,'MinPeakDistance',P1); %'Annotate','extents','WidthReference','halfheight');
[pks8,locs8] = findpeaks(y2,x2,'MinPeakDistance',P1); %'Annotate','extents','WidthReference','halfheight');
hold on
%plot(y2,x2)
plot(locs8,pks8,'x')
legend('Location','northwest');
%Py = meanyr1+stdevyr1; %sets MinPeakProminence value to >stdev
%[pks8,locs8] = findpeaks(y2,x2,'Annotate','extents','MinPeakDistance',Py)
%findpeaks(y1,x2,'Annotate','extents','MinPeakDistance',Py)
%legend('Location','southwest');
xlabel(xlbl1{1});
ylabel(ylbl1{1});
xlim auto;
ylim auto;
locs8mean = mean(diff(locs8))
locs8stdev = std(diff(locs8))
title(['\fontsize{10}Peak-filtered ',num2str(meanyr1x), '-yr ID ',seriesname1{1},' >Mean+2s dist; Mean peak diff = ',num2str(locs8mean),'±',num2str(locs8stdev),' yrs']);

%series2
ax(2) = subplot(2,1,2);
P2 = meanyr2+(2*stdevyr2); %sets MinPeakProminence value to >stdev
findpeaks(y4,x4,'MinPeakDistance',P2); %'Annotate','extents','WidthReference','halfheight');
[pks9,locs9] = findpeaks(y4,x4,'MinPeakDistance',P2); %'Annotate','extents','WidthReference','halfheight');
hold on
%plot(y2,x2)
plot(locs9,pks9,'x')
legend('Location','northwest');
%Py = meanyr1+stdevyr1; %sets MinPeakProminence value to >stdev
%[pks8,locs8] = findpeaks(y2,x2,'Annotate','extents','MinPeakDistance',Py)
%findpeaks(y1,x2,'Annotate','extents','MinPeakDistance',Py)
%legend('Location','southwest');
xlabel(xlbl1{1});
ylabel(ylbl1{1});
xlim auto;
ylim auto;
locs9mean = mean(diff(locs9))
locs9stdev = std(diff(locs9))
title(['\fontsize{10}Peak-filtered ',num2str(meanyr1x), '-yr ID ',seriesname2{1},' >Mean+2s distance; Mean peak diff = ',num2str(locs9mean),'±',num2str(locs9stdev),' yrs']);

%% Peak prominence >1/2 width 
figure (42)
%Peak prominence >1/2 width of standardised, detrended and interpolated peak prominence input data 
ax(1) = subplot(2,1,1);
%series 1
%figure 4b half width - 1s filtering
[pks10,locs10] = findpeaks(y2sln,x2,'Annotate','extents');
findpeaks(y2sln,x2,'Annotate','extents');
legend('Location','northwest');
%P = stdev_all/2; %sets MinPeakProminence value to >stdev
%[a,b] = findpeaks(y2,x2,'MinPeakProminence',P)
%findpeaks(y1,x2,'MinPeakProminence',P)
xlabel(xlbl1{1});
ylabel(ylbl1{1});
xlim auto;
ylim auto;
locs10mean = mean(diff(locs10))
locs10stdev = std(diff(locs10))
title(['\fontsize{10}Peak-filtered ',num2str(meanyr1x), '-yr IDZ data ',seriesname1{1},' Peak width >1/2 height; Mean peak diff = ',num2str(locs10mean),'±',num2str(locs10stdev),' yrs']);

%series2
ax(1) = subplot(2,1,2);
[pks11,locs11] = findpeaks(y4sln,x4,'Annotate','extents');
findpeaks(y4sln,x4,'Annotate','extents');
legend('Location','northwest');
%P = stdev_all/2; %sets MinPeakProminence value to >stdev
%[a,b] = findpeaks(y2,x2,'MinPeakProminence',P)
%findpeaks(y1,x2,'MinPeakProminence',P)
xlabel(xlbl2{1});
ylabel(ylbl2{1});
xlim auto
ylim auto
locs11mean = mean(diff(locs9))
locs11stdev = std(diff(locs11))
title(['\fontsize{10}Peak-filtered ',num2str(meanyr2x), '-yr IDZ ' ,seriesname2{1},' Peak width >1/2 height; Mean peak diff = ',num2str(locs11mean),'±',num2str(locs11stdev),' yrs']);

linkaxes([ax(1),ax(2)],'x')
%% PART 5 - CHANGEPOINT ANALYSIS

% most significant changepoint - for interpolted dataset
% based on largest change in the root-mean-square
% after performing element-wise differentiation to remove any slowly 
% varying trends - then adds  zone added at end

% this runs on detrended series 1 and series 2 data overladi by interpolated factor 1 data
% to use standardized data replace
% y1l > y1ln - data points for series 1 detrended
% y2l > y2ln - data points for series 1 detrended & interpolated
% y3l > y1kn - data points for series 2 detrended
% y4l > y4kn - data points for series 2 detrended & interpolated 

figure(51)
% series 1
ax(1) = subplot(2,1,1);
plot(x1,y1ln,'k.','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',2,'LineWidth',0.2), hold on
plot(x2,y2ln,'r','LineWidth',1);
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
plot(x3,y3ln,'k.','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',2,'LineWidth',0.2), hold on
plot(x4,y4ln,'b','LineWidth',1);
i = findchangepts(diff(y4),'Statistic','rms');
ax = gca;
xp2 = [x4(i) ax.XLim([2 2]) x4(i)];
yp2 = ax.YLim([1 1 2 2]);
patch(datenum(xp2),yp2,[.5 .5 .5],'facealpha',0.1);
title(['\fontsize{10} Most signifcant changepoint (RMS) for ',num2str(factor1),'-yr Interpolation of ' seriesname2{1}]);
xlabel(xlbl2{1});
ylabel(ylbl2{1});
%location of most signifcant changepoint - copy paste x10^4
xp1

%% Most Signifcant CP for RMS and then Fig 33 mean & stdev 
% i.e., where the RMS and then mean and standard deviation of the signal change the most

%Reset mean sample interval to average distance between points 
T1a = mean(intv1)
T2a = mean(intv2)

figure (52)
%series 1
ax(1) = subplot(2,1,1);
%Use this line of code for RMS based CP 
%plot(x1,y1,'ko','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',2,'LineWidth',0.2), hold on
[m,n] = findchangepts(x1,'Statistic','rms')
[pks1,loc1] = findchangepts(y1ln,'Statistic','rms')
mycpplot(y1l,'linear',pks1,loc1,T1a), hold on;

%convert the locations from sample numbers to years BP
pkloc1 = pks1*T1a
pkloc1rms = m*T1a
pkloc1x = (min(x1)+(pks1*T1a))

% Resize x axis to plot in cal years
% **** This only works here because t1=0.667 rounded to 1 *** 
% Need a better fix for this 
% Establish 'XTick' Values
xt = get(gca, 'XTick'); 
% Relabel 'XTick' With 'XTickLabel' Values
set(gca, 'XTick', xt, 'XTickLabel', min(x1)+(xt)*meanyr1);
% **** This only works here because t1=0.667 rounded to 1 *** 
xlim auto
ylim auto
xu = get(gca, 'XTick'); 
%set(gca,'Xtick', xu,'XTickLabel', 18000:500:22500)
%xticks([18000 18500 19000 19500 20000 20500 21000 21500])
xlabel(xlbl1{1});
ylabel(ylbl1{1});

ax(1) = subplot(2,1,2);
%series 2
%[m,n] = findchangepts(y2,'Statistic','rms')
%Use this line of code for RMS based CP 
[pks2,loc2] = findchangepts(y3ln,'Statistic','rms')
mycpplot(y3l,'linear',pks2,loc2,T2a), hold on;
%convert the locations from sample numbers to years BP

%convert the locations from sample numbers to years BP
pkloc2 = pks1*T2a
pkloc2rms = m*T2a
pkloc2x = (min(x1)+(pks1*T2a))

% Resize x axis to plot in cal years
% **** This only works here because t1=0.667 rounded to 1 *** 
% Need a better fix for this 
% Establish 'XTick' Values
xt = get(gca, 'XTick'); 
% Relabel 'XTick' With 'XTickLabel' Values
set(gca, 'XTick', xt, 'XTickLabel', min(x3)+(xt)*meanyr2);
% **** This only works here because t1=0.667 rounded to 1 *** 
xlim auto
ylim auto
%xu = get(gca, 'XTick'); 
%set(gca,'Xtick', xu,'XTickLabel', 18000:500:22500)
%xticks([18000 18500 19000 19500 20000 20500 21000 21500])
xlabel(xlbl2{1});
ylabel(ylbl2{1});
hold off;

%linkaxes([ax(1),ax(2)],'x')

%% Most Signifcant CP for mean & stdev combined
% i.e., where the RMS and then mean and standard deviation of the signal change the most

figure(53)

ax(1) = subplot(2,1,1);
%series 1
[m,n] = findchangepts(x1,'Statistic','rms')
[pks3,loc3] = findchangepts(y1ln,'Statistic','std')
mycpplot(y1l,'linear',pks3,loc3,T1a)

%convert the locations from sample numbers to years BP
pkloc3 = pks3*T1a
pkloc3rms = m*T1a
pkloc3x = (min(x2)+(pks3*T1a))

%convert the locations from sample numbers to years BP
pkloc1 = pks1*T1a
pkloc1rms = m*T1a
pkloc1x = (min(x1)+(pks1*T1a))

% Resize x axis to plot in cal years
% **** This only works here because t1=0.667 rounded to 1 *** 
% Need a better fix for this 
% Establish 'XTick' Values
xt = get(gca, 'XTick'); 
% Relabel 'XTick' With 'XTickLabel' Values
set(gca, 'XTick', xt, 'XTickLabel', min(x1)+(xt)*meanyr1);
% **** This only works here because t1=0.667 rounded to 1 *** 
xlim auto
ylim auto
xu = get(gca, 'XTick'); 
%set(gca,'Xtick', xu,'XTickLabel', 18000:500:22500)
%xticks([18000 18500 19000 19500 20000 20500 21000 21500])
xlabel(xlbl1{1});
ylabel(ylbl1{1});

ax(1) = subplot(2,1,2);
%series 2
[o,p] = findchangepts(x3,'Statistic','rms')
[pks4,loc4] = findchangepts(y3ln,'Statistic','std')
mycpplot(y3l,'linear',pks4,loc4,T2a)

%convert the locations from sample numbers to years BP
pkloc4 = pks4*T2a
pkloc4rms = o*T2a
pkloc4x = (min(x4)+(pks3*T2a))

%convert the locations from sample numbers to years BP
pkloc1 = pks1*T2a
pkloc1rms = m*T2a
pkloc1x = (min(x1)+(pks1*T2a))

% Resize x axis to plot in cal years
% **** This only works here because t1=0.667 rounded to 1 *** 
% Need a better fix for this 
% Establish 'XTick' Values
xt = get(gca, 'XTick'); 
% Relabel 'XTick' With 'XTickLabel' Values
set(gca, 'XTick', xt, 'XTickLabel', min(x3)+(xt)*meanyr2);
% **** This only works here because t1=0.667 rounded to 1 *** 
xlim auto
ylim auto
xu = get(gca, 'XTick'); 
%set(gca,'Xtick', xu,'XTickLabel', 18000:500:22500)
%xticks([18000 18500 19000 19500 20000 20500 21000 21500])
xlabel(xlbl2{1});
ylabel(ylbl2{1});

%% Number of most statistically signifcant changepoints based on RMS (Variance) 
figure (54)
ax(1) = subplot(2,1,1);
%series 1
%numc = chosen number of most signifcant/abrupt changes in variance of the signal
numc = 4
maxnumchg = 10
%findchangepts(v1,'MaxNumChanges',numc)
[r,s] = findchangepts(x1,'Statistic','rms')
[pks5,loc5] = findchangepts(y1l,'Statistic', 'rms', 'MaxNumChanges',maxnumchg)
mycpplot(y1ln,'linear',pks5,loc5,T1a)
pkloc5 = pks5*T1a
% pkloc5rms = r*T1
% equivalent changepoint x-axis location in years 
pkloc5x = (min(x1)+pkloc5)

%convert the locations from sample numbers to years BP
pkloc1 = pks1*T1a
pkloc1rms = m*T1a
pkloc1x = (min(x1)+(pks1*T1a))

% Resize x axis to plot in cal years
% **** This only works here because t1=0.667 rounded to 1 *** 
% Need a better fix for this 
% Establish 'XTick' Values
xt = get(gca, 'XTick'); 
% Relabel 'XTick' With 'XTickLabel' Values
set(gca, 'XTick', xt, 'XTickLabel', min(x1)+(xt)*meanyr1);
% **** This only works here because t1=0.667 rounded to 1 *** 
xlim auto
ylim auto
xu = get(gca, 'XTick'); 
%set(gca,'Xtick', xu,'XTickLabel', 18000:500:22500)
%xticks([18000 18500 19000 19500 20000 20500 21000 21500])
xlabel(xlbl1{1});
ylabel(ylbl1{1});

%series 2
ax(2) = subplot(2,1,2);
%numc = chosen number of most signifcant/abrupt changes in variance of the signal
numc = 4
maxnumchg = 10
%findchangepts(v1,'MaxNumChanges',numc)
[r,s] = findchangepts(x3,'Statistic','rms')
[pks5,loc5] = findchangepts(y3ln,'Statistic', 'rms', 'MaxNumChanges',maxnumchg)
mycpplot(y3l,'linear',pks5,loc5,T2a)
% pkloc5 = pks5*T1
% pkloc5rms = r*T1
% equivalent changepoint x-axis location in years 
% pkloc5x = (min(x3)+pkloc5)

%convert the locations from sample numbers to years BP
pkloc1 = pks1*T2a
pkloc1rms = m*T2a
pkloc1x = (min(x3)+(pks1*T2a))

% Resize x axis to plot in cal years
% **** This only works here because t1=0.667 rounded to 1 *** 
% Need a better fix for this 
% Establish 'XTick' Values
xt = get(gca, 'XTick'); 
% Relabel 'XTick' With 'XTickLabel' Values
set(gca, 'XTick', xt, 'XTickLabel', min(x3)+(xt)*meanyr2);
% **** This only works here because t1=0.667 rounded to 1 *** 
xlim auto
ylim auto
xu = get(gca, 'XTick'); 
%set(gca,'Xtick', xu,'XTickLabel', 18000:500:22500)
%xticks([18000 18500 19000 19500 20000 20500 21000 21500])
xlabel(xlbl2{1});
ylabel(ylbl2{1});

% adds variable to title - removes the residual value which can be useful 
% title([,num2str(numc),' most signifcant/abrupt CPs (Variance)']);

%% Changepoints where the mean and the slope of the signal change most abruptly beyond minimum threshold
% changepoints where the mean and the slope of the signal change most abruptly. 
% Specified minimum improvement is greater then the interpolated mean
%Uses y1 and y3 for series 1 and 2 - which are detrended but not interpolated data

%Series 1
figure (55)
findchangepts(y1ln,'Statistic','linear','MinThreshold',T1a)
[pks6,locs6] = findchangepts(y1ln,'Statistic','linear','MinThreshold',T1a)
%findchangepts(y2,'Statistic','linear', 'MinThreshold', thhold, 'MinDistance',T1/2) %2*stdev_all)
mycpplot(y1ln,'linear',pks6,locs6,T1a)
pkloc6=pks6*T1a
pkloc6x = (min(x1)+pkloc6)

%Resize x axis to plot in cal years
% **** This only works here because t1=0.667 rounded to 1 *** 
% Need a better fix for this 
% Establish 'XTick' Values
xt = get(gca, 'XTick'); 
% Relabel 'XTick' With 'XTickLabel' Values
set(gca, 'XTick', xt, 'XTickLabel', min(x1)+(xt)*meanyr1);
% **** This only works here because t1=0.667 rounded to 1 *** 
xlim auto
ylim auto
xu = get(gca, 'XTick'); 
%set(gca,'Xtick', xu,'XTickLabel', 18000:500:22500)
%xticks([18000 18500 19000 19500 20000 20500 21000 21500])
xlabel(xlbl1{1});
ylabel(ylbl1{1});

%title(['d) Most Signifcant CP (Mean & slope): Min Dist=', num2str(1000*T1)])%adds variable to title

%output data to copy if needed
T1
pks6
pkloc6x

%% Series 2
figure (56)
findchangepts(y3ln,'Statistic','linear','MinThreshold',T2a)
[pks7,locs7] = findchangepts(y3ln,'Statistic','linear','MinThreshold',T2a)
%findchangepts(y2,'Statistic','linear', 'MinThreshold', thhold, 'MinDistance',T1/2) %2*stdev_all)
mycpplot(y3ln,'linear',pks7,locs7,T2a)
pkloc7=pks7*T2a
pkloc7x = (min(x3)+pkloc7)

%Resize x axis to plot in cal years
% **** This only works here because t1=0.667 rounded to 1 *** 
% Need a better fix for this 
% Establish 'XTick' Values
xt = get(gca, 'XTick'); 
% Relabel 'XTick' With 'XTickLabel' Values
set(gca, 'XTick', xt, 'XTickLabel', min(x3)+(xt)*meanyr2);
% **** This only works here because t1=0.667 rounded to 1 *** 
xlim auto
ylim auto
xu = get(gca, 'XTick'); 
%set(gca,'Xtick', xu,'XTickLabel', 18000:500:22500)
%xticks([18000 18500 19000 19500 20000 20500 21000 21500])
xlabel(xlbl1{1});
ylabel(ylbl1{1});

%title(['d) Most Signifcant CP (Mean & slope): Min Dist=', num2str(1000*T1)])%adds variable to title

%output data to copy if needed
T1
pks7 %sample number
pkloc7x %sample number in cal yr BP

%% PART 6 - WAVELET ANALYSIS
%% Example
% This example illustrates how simple it is to do
% continuous wavelet transform (CWT), Cross wavelet transform (XWT)
% and Wavelet Coherence (WTC) plots of your own data.
%
% The time series we will be analyzing are the winter
% Arctic Oscillation index (AO) and
% the maximum sea ice extent in the Baltic (BMI).
%
% http://www.glaciology.net/wavelet-coherence
% close all

%FFT and XWT  Comparison of two data sets,
% should only need to change information at very start of file and plot
% labels
% may also need to alter axis of fourier transforms before saving images
% clear; clc; clf; close all;

%% Choose paths depending on which computer working on
% Add paths for functions and data needed in this sheet on Macbook
% to set cache location
%addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/Wavelet/grinsted-wavelet-coherence-d987ea4');
%addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/Wavelet/grinsted-wavelet-coherence-d987ea4/data');
%addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis');
%addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/PeakAnalysis/Data')
%addpath('/Users/Steve/Dropbox/BAS/Data/MATLAB/Matlab_projects/fLOESS')
% clear; clc; clf; close all; %if needed to clear all 
%% Load the data to run this section independently of previous sections 
% To run wavelet from existing datasets
% First load up the two equally spaced time series into the matrices d1 and d2.

%seriesname={'SAm Temp An' 'Fan-Ti'};
%d1=load('data\SAm.txt');
%d2=load('data\FanTi.txt');

%then goto to section where d1 and d2 are inputs to start 

% or - if not equally spaced time series do below - this is same as before 
% so can skip this section if already run parts 1 to 4

% Import and check 2 data series with age in first column 
% Need to change series loading names and 
% and x and y axis label for plots


% Series 1
%load LP08_I_Mn_Ti.txt
% Series 2
%load LP08_I_Fe_Mn.txt

% Comparing two series
% Create common file name and take natural log of the data
% y = log(x) returns the natural logarithm ln(x) of each element in the array
% change filename to match above

%format long;
%series1x = (LP08_I_Mn_Ti(:,1));
%series1y = log(LP08_I_Mn_Ti(:,2));

% check output - if needed - ; stops the output displaying in command window
%series1x;
%series1y;

% Series 2
%format long;
%series2x = (LP08_I_Fe_Mn(:,1));
%series2y = log(LP08_I_Fe_Mn(:,2));

% Create series labels and x, y axis labels for all plots
seriesname1=seriesname1
%{'LP08 Ln(Mn/Ti)'};
%ylbl1 = {'Ln(Mn/Ti)'};
%xlbl1 = {'Age (cal yr BP)'};

seriesname2=seriesname2
%{'LP08 Ln(Fe/Mn)'};
%ylbl2 = {'Ln(Fe/Mn)'};
%xlbl2 = {'Age (cal yr BP)'};

%% Creating evenly-spaced time series using linear interpolation
% Comparing two time series - need to create even time spacing first
% using detrended y series data 

%series 1 data 
%intv1 = diff(series1x);
%factor1 = 10;
%T1 = round(mean(intv1)*factor1); %dont use this - confusing with precious
%use of T1

% eg here T1 & T2 = 0.667 - therefore rounding = 1 yr x10 = min of 10 yr
% as before, this helps to remove annual-interannual boise in the dataset
% x1 = series1x; %loads x value from column 1, detrends and renames x1
% y1 = detrend(series1y); %loads y value from column2, detrends and renames y1
% Fs1 = 1/T1; %not interpolated  

% Using standardised detrended and PCHIP interpolated data for series 1
xq1 = (yrmin1:meanyr1x:yrmax1)';
vq1 = interp1(x1,y1ln,xq1,'pchip');

%xq1 = x1
%vq1 = y1

%% Series 2 data 
%intv1 = diff(series1x);
%factor2 = 10; 
% eg here T1 & T2 = 0.667 - therefore rounding = 1 yr x10 = min of 10 yr
% as before, this helps to remove annual-interannual boise in the dataset
%T2 = round(mean(intv1)*factor2);
%x2 = series2x; %loads x value from column 1 and renames x1
%y2 = detrend(series2y); %loads y value from column2 and renames y1
%Fs2 = 1/T2; 

% Using standardised detrended and PCHIP interpolated data for series 2
xq2 = (yrmin2:meanyr2x:yrmax2)';
vq2 = interp1(x3,y3ln,xq2,'pchip');

%xq2 = x3
%vq2 = y3

%% Plotting pchip interpolated time series on same figure

figure(61)
%Series 1
ax(1) = subplot(2,1,1);
plot(x1,y1ln,'ko','MarkerSize',1), hold on
plot(xq1,vq1,'r','LineWidth',1), %hold off
xlabel(xlbl1{1});
ylabel(ylbl1{1});
yticks ([-10 -5 -2 -1 0 1 2 5 10]);
title(['\fontsize{12}',num2str(meanyr1x),'-yr IDZ (LBF) ' seriesname1{1}]);

%Series 2
ax(2) = subplot(2,1,2);
plot(x3,y3ln,'ko','MarkerSize',1), hold on
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

%% Change timeseries into percentiles if highly bimodal - not usually
% The time series of Baltic Sea ice extent is highly bi-modal and we
% therefore transform the timeseries into a series of percentiles. The
% transformed series probably reacts 'more linearly' to climate.
%d2(:,2)=boxpdf(d2(:,2));

%% Setting up wavelet plots
seriesname={[seriesname1{1}],[seriesname2{1}]}; %'Ln(Mn/Ti)' 'Ln(Fe/Mn)'};
d1 = [xq1,vq1];
d2 = [xq2,vq2];

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
title(['\fontsize{10}',num2str(meanyr1x),'-yr interpolated Wavelet Power Spectrum of S-D ', seriesname{1}]);
set(gca,'xlim',tlim);

subplot(2,1,2)
wt(d2)
title(['\fontsize{10}',num2str(meanyr2x),'-yr interpolated Wavelet Power Spectrum of S-D ', seriesname{2}])
set(gca,'xlim',tlim)


%% Cross wavelet transform (XWT)
% The XWT finds regions in time frequency space where
% the time series show high common power.
% 'The 5% significance level against red noise is shown as a thick contour. 
% The relative phase relationship is shown as arrows 
% with in-phase pointing right, anti-phase pointing left' 
% where one series leads another by 90?, arrows point straight down).

figure (63) %('color',[1 1 1])
xwt(d1,d2)
title(['\fontsize{10}',num2str(factor1),'-yr interpolated Cross Wavelet Transform (XWT) of S-D ' seriesname{1} ' vs ' seriesname{2} ] )

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
title(['\fontsize{10}',num2str(factor1),'-yr interpolated Wavelet Transform Coherence (WTC) of S-D ' seriesname{1} ' vs ' seriesname{2} ] );

%% Wavelet coherence (WTC) plotted agianst original/interpolated time series
% The WTC finds regions in time frequency space where the two
% time series co-vary (but does not necessarily have high power).

figure(65) %('color',[1 1 1])
ax(1) = subplot(3,1,1);
wtc(d1,d2)
title(['\fontsize{10}',num2str(factor1),'-yr Interpoled Wavelet Transform Coherence (WTC) of Std-Det: ' seriesname{1} ' vs ' seriesname{2} ] );

ax(2) = subplot(3,1,2);
plot(x1,y1,'ko','MarkerSize',0.5), hold on
plot(xq1,vq1,'r','LineWidth',1), %hold off
ylabel(ylbl1{1});

ax(3) = subplot(3,1,3);
plot(x3,y3,'ko','MarkerSize',0.5), hold on
plot(xq2,vq2,'r','LineWidth',1), hold off
xlabel(xlbl2{1});
ylabel(ylbl2{1});

%link parts a to d and zoom in to show the changes
linkaxes(ax(1:3),'x')
xlim auto
ylim auto

%% Copyright notice
%   Copyright (C) 2002-2004, Aslak Grinsted
%
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.


%% This section isn't working - need to debug
figure (26)
%Series 1 - calculate lag time for x-axis of the plot as half the mean interval
lag1 = (meanyr1x)/2;
lagyrs1 = (lag1*50);
lagerr1 = ((std(intv1))/2);
[xc1,lags1] = xcorr(residuals1,lagyrs1,'coeff');%finds autocorrelation with lag in years
%display (xc1);%displays the residuals
%find the critical values at 95% conf level
conf95 = sqrt(2)*erfcinv(0.95);
%conf99 = sqrt(2)*erfcinv(2*.01/2);
lconf1 = -conf95/sqrt(length(x2));
upconf1 = conf95/sqrt(length(x2));
xlim auto;
ylim auto;

% plot the autocorrelation sequence along with the 95%-confidence intervals.
% enter lag time in years to plot in figure 3 in figure 2 where shown

ax(1) = subplot(2,1,1);
stem(lags1,xc1,'filled', 'markersize', 2);
hold on
plot(lags1,[lconf1;upconf1]*ones(size(lags1)),'r','linewidth',2);
ylim auto %ylim([lconf1-0.5 1]);
xlim auto %([-400 400])%auto
xlabel('Lag time (years)');
ylabel('Correlation coefficient');
plot(lags1,lconf1*ones(size(lags1)),'r','linewidth',1);
plot(lags1,upconf1*ones(size(lags1)),'r','linewidth',1);
title(['\fontsize{10} Sample AC with 95% CI for ',seriesname1{1}]);
 
%This doesnt plot same thing as above
%find and plot short and long AC peak locations in residuals
%[pksh1,lcsh1] = findpeaks(xc1);
%short1 = mean(diff(lcsh1))/factor
%[pklg1,lclg1] = findpeaks(xc1,'MinPeakDistance',ceil(short1)*factor*50,'MinPeakHeight',0.1);
%long1 = mean(diff(lclg1))/factor
%pks = plot(lags1(lcsh1)/factor,pksh1,'or', ...
    %lags1(lclg1)/factor,pklg1+0.05,'vk');
%hold off
%legend(pks,[repmat('Period: ',[2 1]) num2str([short2;long2],0)]);
%legend('Location','southwest')
%xlim ([-500 500])%auto
%ylim ([-0.6 1.2])%auto

%% This dosnt work because residuals are = 0 after detrending 

%Series 2 - calculate lag time for x-axis of the plot as half the mean interval
lag2 = (meanyr2x)/2;
lagyrs2 = (lag2*50); %x by 50 to see long term patterns over several cycles 
lagerr2 = ((std(intv2))/2);
[xc2,lags2] = xcorr(residuals2,lagyrs2,'coeff'); %finds autocorrelation with lag in years
% display (xc2); %displays the residuals
%find the critical values at 95% conf level
conf95 = sqrt(2)*erfcinv(0.95);
conf99 = sqrt(2)*erfcinv(2*.01/2);
lconf2 = -conf95/sqrt(length(x4));
upconf2 = conf95/sqrt(length(x4));

ax(2) = subplot(2,1,1);
stem(lags2,xc2,'filled', 'markersize', 2);
hold on
plot(lags2,[lconf2;upconf2]*ones(size(lags2)),'r','linewidth',2);
ylim auto
xlim auto %([-400 400])%auto
xlabel('Lag time (years)');
ylabel('Correlation coefficient');
plot(lags,lconf*ones(size(lags)),'r','linewidth',1);
plot(lags,upconf*ones(size(lags)),'r','linewidth',1);
title(['\fontsize{10} Sample AC with 95% CI for ',seriesname2{1}]);
plot(lags2,[lconf2;upconf2]*ones(size(lags2)),'r','linewidth',2);
hold off

%find and plot peak locations for autocorrelation signal in residuals
%[pksh2,lcsh2] = findpeaks(xc2);
%short2 = mean(diff(lcsh2))/factor;
%[pklg2,lclg2] = findpeaks(xc2,'MinPeakDistance',ceil(short2)*factor*50,'MinPeakHeight',0.1);
%long2 = mean(diff(lclg2))/factor;
%pks = plot(lags2(lcsh2)/factor,pksh2,'or', ...
    %lags2(lclg2)/factor,pklg2+0.05,'vk');
%hold off
%legend(pks,[repmat('Period: ',[2 1]) num2str([short2;long2],0)]);
%legend('Location','southwest')
%xlim ([-50 50])%auto
%ylim ([-0.6 1.2])%auto

% AUTOCORRELATION FIGURE INTERPRETATION

% **** If and when the autocorrelation sequence lies entirely within the 95% CI 
%you can conclude that the residuals are white noise ***

%If the autocorrelation sequence of the residuals looks like the autocorrelation of a white noise process,
%you are confident that none of the signal has escaped your fit and ended up in the residuals.
%Check if and where autocorrelation values clearly exceed the 95%-confidence 
%bounds for a white noise autocorrelation. 
%If they do reject the hypothesis that the residuals are a white noise sequence. 
%The implication is that the model has not accounted for all the signal and 
%therefore the residuals consist of signal plus noise indicating serial
%autocorrelation in the sequence

 