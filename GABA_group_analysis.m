%% GABA-AD: GROUP STATISTICS
% Written by Dominika for GABA-AD project (2021)
% 
% ----- DATA EXTRACTION -----
% Extracts subject info, parameters and outcome measures from all available 
% datasets and saves it in structure 'GABA_results'
% 1) subject info
%       - prepares outcome structure: one row per participant
%       - fills in fields 'subject' and 'session' --> rMT, vigilance etc.
% 2) TEP
%       - dataset info --> number of ICs removed during ICA
%       - mean peak amplitudes and latencies
%       - change in peak amplitude 
%       - SICI
% 3) RS-EEG
%       - dataset info --> ICs, epochs etc.
%       - data parameters --> IAF, TF, fbands 
%       - mean amplitude over fbands
%       - amplitude change
%       - alpha attenuation coeficient
%       - spectral exponent beta 
% 4) MEP
% 
% ----- STATISTICS -----
% 5)



%% 1) prepare outcome structure, subject info
clear all; clc

% parameters
participant = 1:20;
params = readtable('YC_parameters.csv');

% fill in basic info
GABA_results = struct;
for p = 1:length(participant)
    % subject info
    GABA_results(p).subject.age = params{p, 'Age'}; 
    GABA_results(p).subject.sex = params{p, 'Sex'};                         % 1 = female; 0 = male
    GABA_results(p).subject.handedness = params{p, 'Handedness'};           % 1 = right; 0 = left
    GABA_results(p).subject.APOE = params{p, 'APOE'};                       % 1 = homo/hetero eta4; 0 = no eta4
    
    % session info
    GABA_results(p).session(1).medication = 'placebo';
    GABA_results(p).session(1).rmt.pre = params{p, 'zyrtec_rMT_pre'}; 
    GABA_results(p).session(1).rmt.post = params{p, 'zyrtec_rMT_post'}; 
    GABA_results(p).session(1).rmt.change = params{p, 'zyrtec_rMT_post'} - params{p, 'xanax_rMT_pre'}; 
    for a = 1:3 
        if strcmp(class(params{p, 14 + a}), 'cell')
            GABA_results(p).session(1).vigilance(a) = 0;                    % 0 = missing value
        else
            GABA_results(p).session(1).vigilance(a) = params{p, 14 + a};
        end
    end
    
    GABA_results(p).session(2).medication = 'alprazolam';
    GABA_results(p).session(2).rmt.pre = params{p, 'xanax_rMT_pre'}; 
    GABA_results(p).session(2).rmt.post = params{p, 'xanax_rMT_post'}; 
    GABA_results(p).session(2).rmt.change = params{p, 'xanax_rMT_post'} - params{p, 'xanax_rMT_pre'}; 
    for a = 1:3 
        if strcmp(class(params{p, 7 + a}), 'cell')
            GABA_results(p).session(2).vigilance(a) = 0;
        else
            GABA_results(p).session(2).vigilance(a) = params{p, 7 + a};
        end
    end
end

% save
save('GABA_results.mat', 'GABA_results')

%% 2) extract TEP measures
% parameters
stimulus = {'CS' 'TS' 'ppTMS'};
medication = {'placebo' 'alprazolam'};
time = {'pre' 'post'};

% fill in TEP related info
for p = 1:length(participant)
    GABA_results(p).TEP(1).medication = 'placebo';
    GABA_results(p).TEP(1).ICA.pre = params{p, 'zyrtec_ICA_pre'};
    GABA_results(p).TEP(1).ICA.post = params{p, 'zyrtec_ICA_post'};
    
    GABA_results(p).TEP(2).medication = 'alprazolam';
    GABA_results(p).TEP(2).ICA.pre = params{p, 'xanax_ICA_pre'};
    GABA_results(p).TEP(2).ICA.post = params{p, 'xanax_ICA_post'};
end

% mean peak amplitudes and latencies
load('TEPs_final.mat')
for p = 1:length(participant)
    for m = 1:length(medication)
        for s = 1:length(stimulus)
            for t = 1:length(time)
                % choose data --> peak amplitudes and latencies
                rows = (TEPs_final.subject == p & ...
                    categorical(TEPs_final.medication) == medication{m} & ...
                    categorical(TEPs_final.stimulus) == stimulus{s} & ...
                    categorical(TEPs_final.time) == time{t});
                amps = TEPs_final.amplitude(rows)';
                lats = TEPs_final.latency(rows)';
                
                % append to the structure
                statement = ['GABA_results(p).TEP(m).' stimulus{s} '.' time{t} '(1, :) = amps;'];
                eval(statement)
                statement = ['GABA_results(p).TEP(m).' stimulus{s} '.' time{t} '(2, :) = lats;'];
                eval(statement)
            end
        end
    end
end
clear rows amps lats statement 

% amplitude change
TEP_change = readtable('TEP_change_long.csv');
for p = 1:length(participant)
    for m = 1:length(medication)
        for s = 1:length(stimulus)
            for t = 1:length(time)
                % choose data 
                rows = (TEP_change.subject == p & ...
                    categorical(TEP_change.medication) == medication{m} & ...
                    categorical(TEP_change.stimulus) == stimulus{s});
                amp_change = TEP_change.amplitude(rows)';
                
                % append to the structure
                statement = ['GABA_results(p).TEP(m).' stimulus{s} '.change(1, :) = amp_change;'];
                eval(statement)
            end
        end
    end
end
clear rows amp_change statement 

% SICI

% save
save('GABA_results.mat', 'GABA_results')

%% 3) extract RS-EEG measures
% parameters
condition = {'open' 'closed'};

% dataset info --> ICs, epochs etc.

% data parameters --> IAF, TF, fbands 

% mean amplitude over fbands

% amplitude change

% alpha attenuation coeficient

% spectral exponent beta 

% save
save('GABA_results.mat', 'GABA_results')

%% 4) extract MEP measures
















