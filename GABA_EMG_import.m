%% IMPORT EMG DATA TO LETSWAVE
% Written by Dominika for the GABA-AD project (2020)
% 
% 1) Loads raw EMG datasets saved in VHDR format
%    and saves them in the current directory as .mat data + .lw6 header
% 2) Discards 'Out' events, keeps 'Stimulation' and counts them

%% session info
clear all; clc;

% timepoint relative to medication
time = {'pre' 'post'}; 
block = 1:3;

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

        % save the data and the header as letswave files
        header.name = dataset_name;
        save([dataset_name '.mat'], 'data');
        save([dataset_name '.lw6'], 'header');
    end
end