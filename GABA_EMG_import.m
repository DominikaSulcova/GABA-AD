%% IMPORT EMG DATA TO LETSWAVE
% Written by Dominika for the GABA-AD project (2020)
% 
% 1) Loads raw EMG datasets saved in VHDR format
% 2) Keeps only event markers eith the target eventcode (for MEGA --> 's1')
%       - check for duplicate events (identical event latency) and discard
%       repeats
% 3) Saves in the current directory as .mat data + .lw6 header

%% session info
clear all; clc;

% timepoint relative to medication
time = {'pre' 'post'}; 
block = 1:3;
eventcode = 's1';

% subject group
session_info{1} = questdlg('Choose subject group :', 'Aged subjects',...
    'MCI', 'MCI-CTRL', 'CNRAD', 'none');
switch session_info{1}
    case 'MCI'
        group_n = 1;
    case 'MCI-CTRL'
        group_n = 2;
    case 'CNRAD'
        group_n = 3;
end

% subject number
prompt = {'Subject number :'};
dlgtitle = 'Session info';
dims = [1 50];
definput = {[num2str(group_n) '00']};
session_info{2} = char(inputdlg(prompt,dlgtitle,dims,definput));
clear prompt dlgtitle dims definput group_n

% create a prefix
prefix = ['EMG ' session_info{1} ' ' session_info{2}]; 

% choose the folder with raw data
path = sprintf('E:\\Data\\APS-TEP - data\\Raw data\\AGED\\%s%s\\EMG', session_info{1}, session_info{2});
input_folder = uigetdir(path, 'choose folder');
for a = 1:length(time)    
    EMG_code{a} = uigetfile([path '\*.vhdr'], [time{a} ' medication: select file name']);
end
clear a path

%% import MEGA datasets
% import the datasets - blocks indicated by folders vector
for a = 1:length(time)
    for b = 1:legth(block)
        % import the appropriate dataset
        filename = [input_folder '\' time{a} '\b' num2str(block(b)) '\' EMG_code{a}];
        [header, data] = EMG_import_VHDR(filename);

        % create the name for the dataset
        dataset_name = [prefix ' ' time{a} ' b' num2str(block(b)];

        % create the first history entry
        load('EMG_history_import.mat')
        EMG_history_import.configuration.parameters.filenames =  filename;
        header.history(1) =  EMG_history_import;
        
        % keep only the events with target code 
        for c = 1:length(header.events)
            if strcmp(header.events(c).code, eventcode)
                index(c) = true;
            else
                index(c) = false;
            end
        end
        header.events = header.events(index);

        % remove repeated events (based on latency)
        events_unique = unique(extractfield(header.events, 'latency'));
        index_rep = [];
        for d = 1:length(events_unique)
            e = find(extractfield(header.events, 'latency') == events_unique(d));
            index_rep(end + 1) = e(1);
        end
        header.events = header.events(index_rep);

        % save the data and the header as letswave files
        header.name = dataset_name;
        save([dataset_name '.mat'], 'data');
        save([dataset_name '.lw6'], 'header');
    end
end
clear a b c d e events_unique index index_rep