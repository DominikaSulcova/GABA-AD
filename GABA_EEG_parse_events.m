%% PARSE EVENTS ACCORDING TO EXPERIMENTAL CONDITIONS
% Written by Domi for GABA-AD project (2019)
% 
% 1) Loads the order of events for each participant (.xlsx file)
% 2) Loads header for each relevant dataset, assigns codes to events, saves 
% 
% Important:   
%   - Check whether all the files are disponible in the current directory 
%   - Check that all redundent categories of events have been deleted 
%       --> keep 'Stimulation' category for EEG data and 's1' category for EMG data 
%% parameters
clear all
clc

subj_group = 'CNRAD';
participant = [302 303 306 307 308 310];
time = {'pre' 'post'};
block = 1:3;

prefix_EEG = 'reref ds art-sup ep dc YC';
prefix_EMG = 'bl ep notch but dc YC';

%% parse events in all relevant datafiles
% switch according to the subject group
switch subj_group
    case 'MCI'
        sg = 'P';
    case 'MCI-CTRL'
        sg = 'C';
    case 'CNRAD'
        sg = 'R';
end

% loop through participants
for p=1:length(participant)
    % create a conditions matrix 
    conditions = xlsread(['conditions_order_' sg '_' num2str(participant(p)) '.xlsx']);
    
    % loop through datasets
    for t=1:length(time)
        % parse EEG events
        for b=1:length(block)
            load([prefix_EEG num2str(participant(p)) ' ' char(session(s)) ' ' char(time(t)) ' M1 block' num2str(block(b)) '.lw6'],'-mat'); 
            for e=1:length(header.events)
                header.events(e).code = num2str(conditions(e, b));
            end
            save([prefix_EEG num2str(participant(p)) ' ' char(session(s)) ' ' char(time(t)) ' M1 block' num2str(block(b)) '.lw6'],'header');
        end

        % parse EMG events
        for b=1:length(block)
            load([prefix_EMG num2str(participant(p)) ' ' char(session(s)) ' ' char(time(t)) ' EMG block' num2str(block(b)) '.lw6'],'-mat'); 
            for e=1:length(header.events)
                header.events(e).code = num2str(conditions(e, b));
            end
            save([prefix_EMG num2str(participant(p)) ' ' char(session(s)) ' ' char(time(t)) ' EMG block' num2str(block(b)) '.lw6'],'header');
        end
        
    end
end