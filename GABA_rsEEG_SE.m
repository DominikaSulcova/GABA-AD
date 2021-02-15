%% RS-EEG - COMPUTE SPECTRAL EXPONENT BETA 
% Written by Dominika for GABA-AD project (2021)
% Based on scripts and functions written by Michele A. Colombo from Massimini's lab, Milano
% 
% 1) Prepares data
%       - Loads FFT transformed data (relative amplitude) - 'rsEEG_data_high.mat'
%           --> medication * time * condition * participant * ROI * datapoint
%       - Log-log transforms of the data

%% parameters
clear all 
clc

% dataset
xstep = 0.25;

% process
window = [1 20; 20 40; 1 40]; 
reg_method = 'ols';

% a = 1; 

%% 1) prepare data, log-log transform
% load the data
load('rsEEG_data_high.mat');

for a = 1:size(window, 1)
    % choose data
    x_start = ; 
    x_end = ;
    x = window(a, 1):xstep:window(a, 2); 
    data = data_high(:, :, :, :, :, x_start : x_end);

    % avoid zeros 
    zero = data == 0; data(zero)=[]; x(zero)=[]; 
end
 
% log transform
x = log10(x);
data_visual_open = log10(data_visual_open);
data_visual_closed = log10(data_visual_closed);

% interpolate in log-log: X, Y --> Xi, Yi
 stepsFrex= (diff( XX)); stepsFrex= round(stepsFrex.*1000)./1000;
 if length( unique(stepsFrex))==1 %IF frex is a vector of linearly increasing frex, try
     XXi= logspace( X(1), X(end),  length(X)*4); % 2^10*2;--> auto-chosen
      Xi= log10(XXi); 
     Yi= interpn( X, Y,  (Xi) );% interpn ,'spline', 'nearest' 
     YYi= 10.^(Yi); 
 else % catch
     XXi= XX;  Xi= log10(XXi); 
     YYi= YY;  Yi= log10(YYi); disp('could not upsample PSD')
 end