%% RS-EEG PREPROCESSING 2
% Written by Dominika for GABA-AD project (2021)
% 
% 1) Reads 'treatments.xlsx' and creates a cell array with medication
%    information for each dataset
% 2) Replaces session information with medication information, saves
% 3) After averaging (done in letswave), merges all subjects in group dataset 
%    for each category --> calls a lw6 function 'RLW_merge_epochs' 

%% parameters
clear all
clc

% dataset
participant = 1:20;  
medication = {'placebo' 'alprazolam'};  
session = {'S1' 'S2'};      
time = {'pre' 'post'};
condition = {'open' 'closed'};
frequency = {'low' 'high'};
prefix_1 = 'han icfilt assigned ica'; 
prefix_2 = 'crop_120 45Hz_but4 but reref chan-select ds YC';

% load and adjust the medication info
med_session = table2cell(readtable('treatments.xlsx'));
for a = 1 : size(med_session, 1)
    for b = 1 : length(session)
        if med_session{a, b+1}(1) == 'z'
            med_session{a, b+1} = 'placebo'; 
        else
            med_session{a, b+1} = 'alprazolam'; 
        end
    end
end

%% replace the session number in the name by corresponding medication
% loop through datasets
% p = 1; s = 1; t = 1; c = 1; f = 1;   
for p = 1:length(participant)
    for s = 1:length(session)
        this_session = med_session(participant(p), s+1);
        for t = 1:length(time)       
            for c = 1:length(condition)
                for f = 1:length(frequency)
                    % rename fully filtered data (all artifactual ICAs removed)
                    name_data = ['fft ' prefix_1 ' ' frequency{f} ' ' prefix_2 num2str(participant(p)) ' ' session{s} ' ' time{t} ' EEGcont ' condition{c} '.mat'];
                    name_header = ['fft ' prefix_1 ' ' frequency{f} ' ' prefix_2 num2str(participant(p)) ' ' session{s} ' ' time{t} ' EEGcont ' condition{c} '.lw6'];
                    load(name_data);
                    load(name_header, '-mat');
                    
                    % save under new name
                    name_new = ['fft ' prefix_1 ' ' frequency{f} ' ' prefix_2 num2str(participant(p)) ' ' char(this_session) ' ' time{t} ' EEGcont ' condition{c}];
                    header.name = name_new;
                    save([name_new '.mat'], 'data');
                    save([name_new '.lw6'], 'header');                        
                end 
            end
        end
    end
    clear name_data name_header name_new data header
    message = ['Participant n.' num2str(participant(p)) ' finished.'];
    disp(message)
end

%% merge epochs of each category
% m = 1; t = 1; c = 1; s = 1; p = 1; 

% create 'merge_idx' variable for epoch merging
merge_idx = 1:numel(participant);

% loop through conditions
for m = 1:numel(medication)
    for t = 1:numel(time)        
        for c = 1:numel(condition)
            disp(['>>>>> Processing : ' medication{m} ' - ' time{t} ' - ' condition{c} ' <<<<<'])
            for f = 1:numel(frequency)
                disp(['..... ' frequency{f} ' frequency .....'])              
                
                % create 'datasets' variable for epoch merging
                name_old = ['avg fft ' prefix_1 ' ' frequency{f} ' ' prefix_2 num2str(participant(p)) ' ' medication{m} ' ' time{t} ' EEGcont ' condition{c}];
                datasets = struct;
                for p = 1:numel(participant)
                    load([name_old '.lw6'], '-mat');
                    datasets(p).header = header;
                    load([name_old '.mat']);
                    datasets(p).data = data;
                end
                clear data header

                % merge datasets
                [header,data,message_string] = RLW_merge_epochs(datasets,merge_idx);
                
                % create new name
                name_new = ['merged fft ' frequency{f} ' ' medication{m} ' ' time{t} ' EEGcont ' condition{c}];

                % save the data
                save([name_new '.mat'], 'data')

                % save the header
                header.name = name_new;
                save([name_new '.lw6'], 'header')
                
                disp('... done.')            
            end
        end
    end    
end
clear data header 
disp('>>>>>>> All done. <<<<<<<')
