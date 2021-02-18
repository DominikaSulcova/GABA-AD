%% RS-EEG - COMPUTE SPECTRAL EXPONENT BETA 
% Written by Dominika for GABA-AD project (2021)
% Based on scripts and functions written by Michele A. Colombo from Massimini's lab, Milano
% 
% 1) Prepares data
%       - averages signal across ROIs
%       - cuts x and data according to target fband windows and saves in a
%         structure --> process.x; process.data
% 2) Fits the power using Michele's function 
%       - in 3 steps
%       --> outcome variables:  intslo(1)
%                               intslo(2)                   
%                               


%% parameters
clear all 
clc

% dataset
load('rsEEG_data_high.mat');
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'}; 
condition = {'open' 'closed'};

% process
do_plot = 1;
xstep = 0.25; 
xoffset = 0.1;
load('colours2.mat')

% outcome structure
spect_exp = struct;
spect_exp(1).band = 'low'; spect_exp(2).band = 'high'; spect_exp(3).band = 'broad'; 
spect_exp(1).window = [1 20]; spect_exp(2).window = [20 40]; spect_exp(3).window = [1 40]; 
spect_exp(1).method = 'ols'; spect_exp(2).method = 'ols'; spect_exp(3).method = 'ols';

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
for a = 1:numel(spect_exp)
    % choose original FFT data: x, data
    x_start = ceil((spect_exp(a).window(1) - xoffset) / xstep); 
    x_end = ceil((spect_exp(a).window(2) - xoffset) / xstep);
    x = spect_exp(a).window(1):xstep:spect_exp(a).window(2); 
    spect_exp(a).x = x;
    spect_exp(a).data = data(:, :, :, :, x_start : x_end);   
end
clear x_start x_end x data
clear a m t c p r i
   
%% 2) fit the power law 
% loop through target fband windows
for a = 1:numel(spect_exp)   
    % choose vector with frequency bins
    x = spect_exp(a).x;
    
    % loop through datasets
    for m = 1:size(spect_exp(a).data, 1)
        for t = 1:size(spect_exp(a).data, 2)
            for c = 1:size(spect_exp(a).data, 3)
                for p = 1:size(spect_exp(a).data, 4)
                    % call data
%                     dataset_name = ['YC' num2str(p) '_' medication{m} '_' time{t} '_' condition{c}];
                    data = squeeze(spect_exp(a).data(m, t, c, p, :))';
                    
                    % choose plot colour
                    col = colours2((m - 1)*2 + t, :);
                    
                    % fit current dataset
                    [intslo, stats, amps] = fitPowerLaw3steps(x, data, spect_exp(a).method, do_plot, col); 
%                     hold on
                    
                    % fill in the outcome structure
                    spect_exp(a).result.intercept(m, t, c, p) = intslo(1);
                    spect_exp(a).result.slope(m, t, c, p) = intslo(2);                    
                    
                    % clear figure
                    clf
                end 
            end
        end
    end
end
clear dataset_name x data col intslo stats amps fig
clear a m t c p i

% save outcome structure
save('xpect_exp.mat', 'spect_exp')

