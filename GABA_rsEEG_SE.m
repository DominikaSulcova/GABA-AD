%% RS-EEG - COMPUTE SPECTRAL EXPONENT BETA 
% Written by Dominika for GABA-AD project (2021)
% Based on scripts and functions written by Michele A. Colombo from Massimini's lab, Milano
% 
% 1) Prepares data
%       - averages signal across ROIs


%% parameters
clear all 
clc

% dataset
load('rsEEG_data_high.mat');
xstep = 0.25; xoffset = 0.1;
medication = {'placeb' 'alprazolam'}; time = {'pre' 'post'}; condition = {'open' 'closed'};
load('colours2.mat')

% process
process = struct;
process(1).band = 'low'; process(2).band = 'high'; process(3).band = 'broad'; 
process(1).window = [1 20]; process(2).window = [20 40]; process(3).window = [1 40]; 
process(1).method = 'ols'; process(2).method = 'ols'; process(3).method = 'ols';

% a = 1; m = 1; t = 1; c = 1; p = 1; r = 1; 

%% 1) prepare data
% average data across regions
for m = 1:size(data_high, 1)
    for t = 1:size(data_high, 2)
        for c = 1:size(data_high, 3)
            for p = 1:size(data_high, 4)
                for i = 1:size(data_high, 6)
                   data(m, t, c, p, i) = mean(data_high(m, t, c, p, :, i));
                end
            end
        end
    end
end
clear data_high

% chops datasets by target fbands
for a = 1:numel(process)
    % choose original FFT data: x, data
    x_start = ceil((process(a).window(1) - xoffset) / xstep); 
    x_end = ceil((process(a).window(2) - xoffset) / xstep);
    x = process(a).window(1):xstep:process(a).window(2); 
    process(a).x = x;
    process(a).data = data(:, :, :, :, x_start : x_end);   
end
clear x_start x_end x
   
%% 2) fit  power law in 3 steps
% loop 
[insl, st, pp] = fitPowerLaw3steps(frex( myX),  myY( myX),[], 0); hold on

