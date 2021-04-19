%% params
clear all, clc

group = {'CNRAD' 'MCI-CTRL' 'MCI'};
participant = {[301 302 303 306 307 308 310] [206] [104]};
time = {'pre' 'post'};
block = 1:3;
prefix = 'EEG M1';
codename = 'Stimulation';

%% rename
for g = 1:length(group)
    for p = 1:length(participant{g})
        for t = 1:length(time)
            for b = 1:length(block)
                load([prefix ' ' group{g} ' ' num2str(participant{g}(p)) ' ' time{t} ' b' num2str(block(b)) '.lw6'], '-mat')
                for e = 1:length(header.events)
                    header.events(e).code = codename;
                end
                save([prefix ' ' group{g} ' ' num2str(participant{g}(p)) ' ' time{t} ' b' num2str(block(b)) '.lw6'], 'header')
            end
        end
    end
end
clear g p t b e

%% concatenate
% params
subject = {'01' '02' '03' '04'};
intensity = {'stim_100' 'stim_120' 'stim_140'};
position = {'across' 'along'};
current = {'normal' 'reversed'}; 
prefix = 'avg avgchan bl icfilt ica visual crop but fft-notchfilt prefilt prea P1';

% loop 
for i = 1:length(intensity)
    for p = 1:length(position)
        for c = 1:length(current)
            % merge datasets
            merge_idx = 1:length(subject);
            datasets = struct;
            for s = 1:length(subject)
                name_old = [prefix ' ' subject{s} ' ' position{p} ' ' current{c} ' ' intensity{i}]; 
                load([name_old '.lw6'], '-mat');
                datasets(s).header = header;
                load([name_old '.mat']);
                datasets(s).data = data;
            end
            [header,data,message_string] = RLW_merge_epochs(datasets,merge_idx); 

            % save datasets
            name_new = ['merged ' position{p} ' ' current{c} ' ' intensity{i}];
            save([name_new '.mat'], 'data')
            header.name = name_new;
            save([name_new '.lw6'], 'header')   
        end
    end
end
clear i c p s merge_idx data header message_string name_new name_old



