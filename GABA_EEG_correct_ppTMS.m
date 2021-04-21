%% ppTMS CORRECTION 
% Written by Dominika for GABA-AD project (2019)
% 
% Corrects TEP elicited by the paired-pulse TMS stimulation for preceding
% CS TEP --> assumes linar summation of evoked potentials, which is ultimately 
% totally faulty but approximately very useful
%   - shifts the average CS dataset backwards by 2.5ms 
%   - substracts the shifted CS dataset from the corresponding ppTMS dataset

%% parameters
clear all
clc

participant = 41:20;
session = {'S1' 'S2'};
time = {'pre' 'post'};
prefix = 'avg avgchan icfilt assigned2 ica26 but fft-notchfilt prefilt assigned prea28 visual complete A10 reref';
operation = 'A-B';

%% core loop
for p=1:length(participant)
    for s=1:length(session)
        for t=1:length(time)
            % deal with the CS data 
            load([prefix ' CS 80 YC' num2str(participant(p)) ' ' char(session(s)) ' ' char(time(t)) ' M1.mat']);
            data = data(:, :, :, :, :, 6:6000);                 % cut away 5 datapoints (2 datapoints = 1 ms) in the beginning
            ending = data(:, :, :, :, :, 1:5);                  % and add 5*zero to the end
            ending(:, :, :, :, :, :) = 0;
            data = cat(6, data, ending);
            save(['shifted ' prefix ' CS 80 YC' num2str(participant(p)) ' ' char(session(s)) ' ' char(time(t)) ' M1.mat'], 'data');
            
            % deal with the header
            load([prefix ' CS 80 YC' num2str(participant(p)) ' ' char(session(s)) ' ' char(time(t)) ' M1.lw6'], '-mat');
            save(['shifted ' prefix ' CS 80 YC' num2str(participant(p)) ' ' char(session(s)) ' ' char(time(t)) ' M1.lw6'], 'header');
            
            % load the ppTMS dataset and perform the substraction
            dataB = data; headerB = header;                     % CS will be dataset B
            load([prefix ' ppTMS 80120 YC' num2str(participant(p)) ' ' char(session(s)) ' ' char(time(t)) ' M1.mat']);
            dataA = data;                                       % ppTMS will be dataset A
            load([prefix ' ppTMS 80120 YC' num2str(participant(p)) ' ' char(session(s)) ' ' char(time(t)) ' M1.lw6'], '-mat');
            headerA = header;
            clear header data ending 
            
            % use a LW function to perform the operation A-B            
            [header,data,message]=RLW_math(headerA,dataA,headerB,dataB,operation);
            
            % save data and header under new name
            header.name = ['corrected ' prefix ' ppTMS 80120 YC' num2str(participant(p)) ' ' char(session(s)) ' ' char(time(t)) ' M1'];
            save([header.name '.mat'], 'data');
            save([header.name '.lw6'], 'header');
            clear header headerA headerB data dataA dataB message
        end
    end 
end