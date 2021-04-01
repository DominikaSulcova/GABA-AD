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
