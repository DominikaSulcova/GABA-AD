%% concatenate
clear all
clc

% params
group = 'YC';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
stimulus = {'CS' 'TS' 'ppTMS'};
time = {'pre' 'post'};
target = 'M1';
prefix = 'eois crop avg final_dataset';
       
% loop 
for m = 1:length(medication)
    for s = 1:length(stimulus)
        for t = 1:length(time)
            % merge datasets
            merge_idx = 1:length(participant);
            datasets = struct;
            for p = 1:length(participant)
                name_old = [prefix ' ' stimulus{s} ' ' group num2str(participant(p)) ' ' medication{m} ' ' time{t} ' ' target]; 
                load([name_old '.lw6'], '-mat');
                datasets(p).header = header;
                load([name_old '.mat']);
                datasets(p).data = data;
            end
            [header,data,message_string] = RLW_merge_epochs(datasets,merge_idx); 

            % save datasets
            name_new = ['merged ' medication{m} ' ' time{t} ' ' stimulus{s} ' ' target];
            save([name_new '.mat'], 'data')
            header.name = name_new;
            save([name_new '.lw6'], 'header')   
        end
    end
end
clear m s t p merge_idx data header message_string name_new name_old

%% rename treatments
clear all
clc

% parameters
participant = 1:20;
stimulus = {'CS 80' 'TS 120' 'ppTMS 80120'};
session = {'S1' 'S2'};

% deal with treatments
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
clear a b 

%% calculate number of epochs
clear all
clc

% parameters
participant = 1:20;
stimulus = {'CS 80' 'TS 120' 'ppTMS 80120'};
session = {'S1' 'S2'};
time = {'pre' 'post'};
prefix = 'long avgchan icfilt assigned2 ica26 but fft-notchfilt prefilt assigned prea28 visual complete A10 reref'; 
path = 'F:\GABA-AD\Data\Processed data\YC\M1_TMS-EEG';

% p = 1; se = 1; t = 1; s = 1;
for p = 1:length(participant)
    for se = 1:length(session)
        % decide session
        this_session = med_session(participant(p), se+1);              
        for t = 1:length(time)
            for s = 1:length(stimulus)
                % load data
                load([path '\' prefix ' ' stimulus{s} ' YC' num2str(participant(p)) ' ' session{se} ' ' time{t} ' M1.mat'])
                
                % extract number of epochs
                if strcmp(this_session{1}, 'placebo')
                    epochs(p, 1, t, s) = size(data, 1);
                elseif strcmp(this_session{1}, 'alprazolam')
                    epochs(p, 2, t, s) = size(data, 1);
                end                
            end
        end
    end
end
clear p se s t data

save('epochs.mat', 'epochs')

%% rename event code
clear all
clc

% parameters
group = 'MCI-CTRL';
participant = [206];
% group = 'CNRAD';
% participant = [302, 306, 307, 308, 310];
time = {'pre' 'post'};
prefix = 'EEG AG'; 
path = 'E:\Data\GABA-AD - data\Processed data\AGED\AG_TMS-EEG';
eventcode = 'Stimulation';

% p = 1; t = 1; e = 1;
for p = 1:length(participant)
    for t = 1:length(time)
        % load header
        load([path '\' prefix ' ' group ' ' num2str(participant(p)) ' ' time{t} '.lw6'], '-mat') 
        
        % rename all events
        for e = 1:length(header.events)
            header.events(e).code = eventcode;
        end
        
        % save the header
        save([path '\' prefix ' ' group ' ' num2str(participant(p)) ' ' time{t} '.lw6'], 'header')
    end
end
clear p t e 


