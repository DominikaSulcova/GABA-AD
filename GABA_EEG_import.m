%% Script to import EEG datasets 
% Author: Dominika
% Study: GABA-AD
% Data: 
% 1) RS-EEG
%       - 3 mins eyes open
%       - 3 mins eyes closed
% 2) TMS stimulation over left M1 
%       - 60 stims 80 %rMT   --> CS 
%       - 60 stims 120 %rMT  --> TS
%       - 60 stims paired pulse (CS - 2.5 ms - TS) --> ppTMS
% 3) TMS stimulation over left AG 
%       - 80 stims 120 %rMT 

% BEFORE THE SCRIPT RUNS:
% - Check specific parameters in the beginning of each section
% - The script mports each session separately and asks for numbers of folders that
% contain the data --> verify 
% - In case the data were registered into several folders per timepoint
% (e.g. when the session is interupted), it is possible to specify new folder 
% at the beginning of each section. Comment the section entitled 'WHEN
% RUNNING THE WHOLE SCRIPT' and enable sections entitled 'WHEN RUNNING
% SECTION BY SECTION'. Then run the script by sections,  obviously :P

%% session info
clear all; clc;

% timepoint relative to medication
time = {'pre' 'post'}; 

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
dims = [1 35];
definput = {[num2str(group_n) '00']};
session_info{2} = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput group_n

%% create a logfile
% get the date
prompt = {'Date of recording :'};
dlgtitle = 'Date';
dims = [1 40];
definput = {'00/00/2020'};
session_info{3} = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput 

% inicialize the logfile
filename = GABA_initialize_logfile(session_info);

%% RS-EEG recording
% dataset parameters
eyes = {'open' 'closed'};
prefix = ['RS-EEG ' char(session_info{1}) ' ' char(session_info{2})];

% choose folders with raw data
folder_path = ['C:\Users\sulcova\Desktop\GABA-AD\' char(session_info{1})];
session_info{4} = uigetdir(folder_path, 'Folder: pre-medication');
session_info{5} = uigetdir(folder_path, 'Folder: post-medication');

% specify EEG block subfolders in the imput folder
prompt = {'Pre-medication' 'Post-medication'};
dlgtitle = 'Admit following EEG blocks :';
dims = [1 60];
definput = {'[1 2]' '[1 2]'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
session_info{6} = cell2num(answer(1));
session_info{7} = cell2num(answer(2));
clear answer prompt dlgtitle dims definput 

% -------------------- WHEN RUNNING THE WHOLE SCRIPT -------------------- 
f = 0;
f2 = 0;

% specify EEG block subfolders in the imput folder
prompt = {'Pre-medication' 'Post-medication'};
dlgtitle = 'Admit following EEG blocks :';
dims = [1 60];
definput = {'[3:5]' '[3:5]'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
session_info{8 + f} = cell2num(answer(1));
session_info{9 + f} = cell2num(answer(2));
clear answer prompt dlgtitle dims definput 

% specify EEG block subfolders in the imput folder
prompt = {'Pre-medication' 'Post-medication'};
dlgtitle = 'Admit following EEG blocks :';
dims = [1 60];
definput = {'6' '6'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
session_info{10 + f + f2} = cell2num(answer(1));
session_info{11 + f + f2} = cell2num(answer(2));
clear answer prompt dlgtitle dims definput 
% -----------------------------------------------------------------------

% loop through timepoints and blocks
for a = 1:length(time)
    for b = 1:length(eyes)
        % import the appropriate dataset
        [header, data] = RLW_import_MEGA(session_info{4 + a - 1}, session_info{6 + a - 1}(b));

        % create the name for the dataset
        dataset_name = [prefix ' ' time{a}  ' ' eyes{b}];

        % create the first letswave history entry
        load('EEG_history_import.mat')
        EEG_history_import.configuration.parameters.input_folder  = session_info{4 + a - 1};
        EEG_history_import.configuration.parameters.session_number  = session_info{6 + a - 1}(b);   
        header.history(1) =  EEG_history_import;

        % save the data and the header as letswave files
        header.name = dataset_name;
        save([dataset_name '.mat'], 'data');
        save([dataset_name '.lw6'], 'header');
    end
end

%% M1 recording
% dataset parameters
block = {'b1' 'b2' 'b3'};
prefix = ['EEG M1 ' char(session_info{1}) ' ' char(session_info{2})];
event_code = 'B - Stimulation';

% % ------------------- WHEN RUNNING SECTION BY SECTION -------------------
% % check for validity of raw data location 
% answer = questdlg('Continue importing from previous data folders?', 'Raw data folders',...
%     'Yes', 'No, change folder location', 'Yes');
% switch answer
%     case 'Yes'
%         f = 0;
%     case 'No, change folder location'
%         f = 2;
%         % choose folders with raw data
%         folder_path = ['C:\Users\sulcova\Desktop\GABA-AD\' char(session_info{1})];
%         session_info{8} = uigetdir(folder_path, 'Folder: pre-medication');
%         session_info{9} = uigetdir(folder_path, 'Folder: post-medication');        
% end
% clear answer
% 
% % specify EEG block subfolders in the imput folder
% prompt = {'Pre-medication' 'Post-medication'};
% dlgtitle = 'Admit following EEG blocks :';
% dims = [1 60];
% definput = {'[3:5]' '[3:5]'};
% answer = inputdlg(prompt,dlgtitle,dims,definput);
% session_info{8 + f} = cell2num(answer(1));
% session_info{9 + f} = cell2num(answer(2));
% clear answer prompt dlgtitle dims definput
% % -----------------------------------------------------------------------

% loop through timepoints and blocks
for a = 1:length(time)
    for b = 1:length(block)
        % import the appropriate dataset
        [header, data] = RLW_import_MEGA(session_info{4 + f * 2 + a - 1}, session_info{8 + f + a - 1}(b));

        % ditch extra 'out' category
        T = struct2table(header.events);                                    % creates a table of the structure 'events'
                    T.code = categorical(T.code);                           % in order to work with strings --> convert to categorical
                    sortedT = T(T.code == event_code, :); 
                    sortedT.code = cellstr(sortedT.code);                   % turns the categorical data back to string cells
                    header.events = table2struct(sortedT);                  % turns the sorted table back in a structure field 
                    header.events = header.events';

        % create the name for the dataset
        dataset_name = [prefix ' ' time{a}  ' ' block{b}];

        % create the first letswave history entry
        load('EEG_history_import.mat')
        EEG_history_import.configuration.parameters.input_folder  = session_info{4 + f * 2 + a - 1};
        EEG_history_import.configuration.parameters.session_number  = session_info{8 + f + a - 1}(b);   
        header.history(1) =  EEG_history_import;

        % save the data and the header as letswave files
        header.name = dataset_name;
        save([dataset_name '.mat'], 'data');
        save([dataset_name '.lw6'], 'header');
    end
end

%% AG recording
% dataset parameters
prefix = ['EEG AG ' char(session_info{1}) ' ' char(session_info{2})];

% ------------------- WHEN RUNNING SECTION BY SECTION -------------------
% % check for validity of raw data location 
% answer = questdlg('Continue importing from previous data folders?', 'Raw data folders',...
%     'Yes', 'No, change folder location', 'Yes');
% switch answer
%     case 'Yes'
%         f2 = 0;
%     case 'No, change folder location'
%         f2 = 2;
%         % choose folders with raw data
%         folder_path = ['C:\Users\sulcova\Desktop\GABA-AD\' char(session_info{1})];
%         session_info{10 + f} = uigetdir(folder_path, 'Folder: pre-medication');
%         session_info{11 + f} = uigetdir(folder_path, 'Folder: post-medication');        
% end
% clear answer

% % specify EEG block subfolders in the imput folder
% prompt = {'Pre-medication' 'Post-medication'};
% dlgtitle = 'Admit following EEG blocks :';
% dims = [1 60];
% definput = {'6' '6'};
% answer = inputdlg(prompt,dlgtitle,dims,definput);
% session_info{10 + f + f2} = cell2num(answer(1));
% session_info{11 + f + f2} = cell2num(answer(2));
% clear answer prompt dlgtitle dims definput 
% -----------------------------------------------------------------------

% loop through timepoints and blocks
for a = 1:length(time)
    % import the appropriate dataset
    [header, data] = RLW_import_MEGA(session_info{4 + f * 2 + f2 * 3 + a - 1}, session_info{10 + f + f2 + a - 1});

    % ditch extra 'out' category
    T = struct2table(header.events);                                    % creates a table of the structure 'events'
                T.code = categorical(T.code);                           % in order to work with strings --> convert to categorical
                sortedT = T(T.code == 'B - Stimulation', :); 
                sortedT.code = cellstr(sortedT.code);                   % turns the categorical data back to string cells
                header.events = table2struct(sortedT);                  % turns the sorted table back in a structure field 
                header.events = header.events';

    % create the name for the dataset
    dataset_name = [prefix ' ' time{a}];

    % create the first letswave history entry
    load('EEG_history_import.mat')
    EEG_history_import.configuration.parameters.input_folder  = session_info{4 + f * 2 + f2 * 3 + a - 1};
    EEG_history_import.configuration.parameters.session_number  = session_info{10 + f + f2 + a - 1};   
    header.history(1) =  EEG_history_import;

    % save the data and the header as letswave files
    header.name = dataset_name;
    save([dataset_name '.mat'], 'data');
    save([dataset_name '.lw6'], 'header');
end