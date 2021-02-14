%% RS-EEG PREPROCESSING 1
% Written by Dominika for GABA-AD project (2021)
% 
% 1) Checks the number of available datasets, in case of split datasets
%    saves the information into a cell array 'split_datasets'
% 2) Resets the first timepoint value to 0 
% 3) After segmentation (performed in letswave), the second section of the
%    script performs merging of possible split datasets for each segmentation category

%% parameters
clear all
clc

participant = 1:20;
session = {'S1' 'S2'};
time = {'pre' 'post'};
condition = {'open' 'closed'};
prefix = '45Hz_but4 but reref chan-select ds YC';
seg_category = {'ep_ica' 'low' 'high'};
split_suffix = {'A' 'B' 'C' 'D' 'E'};

%% reset the timepoint count from 0
% prepare a cell array of split datasets for future merging
split_datasets = {};

% loop through datasets
% p = 1; s = 1; t = 1; c = 1; 
for p = 1:length(participant)
    for s = 1:length(session)
        for t = 1:length(time)
            for c = 1:length(condition)
                % search for available datasets
                base_name = [prefix num2str(participant(p)) ' ' char(session(s)) ' ' char(time(t)) ...
                    ' EEGcont ' char(condition(c)) '.lw6'];
                datasets = dir(['*crop_120*' base_name]);
                
                % if the dataset is split, save basename for future merging of epochs
                if numel(datasets) > 1
                    split_datasets = [split_datasets, base_name];
                end                        
                                
                % looop through datasets
                for d = 1:numel(datasets)
                    % load the header
                    load(datasets(d).name, '-mat')
                    
                    % reset the time count
                    header.xstart = 0;
                    
                    % save the header
                    save(datasets(d).name, 'header');                    
                end 
                disp(['>>> Participant n.' num2str(participant(p)) ' - ' session{s} ' - ' time{t} ' medication - eyes ' condition{c} ' : DONE. <<<'])
            end 
        end
    end
end

%% merge epochs of split datasets
% to be done AFTER segmentation to chunks

% a = 1; b = 1; c = 1; 
for a = 1:numel(split_datasets)
    % define the base name
    base_name = split_datasets{a};
    disp(['>>>>>>> crop_120 ' base_name(1:end-4) ' <<<<<<<'])
        
    % loop through segmentation categories
    for b = 1:numel(seg_category)
        % search for available datasets
        datasets2merge = dir(['*' seg_category{b} '*crop_120*' base_name]);
    
        % create 'merge_idx' variable for epoch merging
        merge_idx = 1:numel(datasets2merge);
        
        % create 'datasets' variable for epoch merging
        datasets = struct;
        for c = 1:numel(datasets2merge)
            load([seg_category{b} ' crop_120_' split_suffix{c} ' ' base_name], '-mat');
            datasets(c).header = header;
            load([seg_category{b} ' crop_120_' split_suffix{c} ' ' base_name(1:end-4) '.mat']);
            datasets(c).data = data;
        end
        
        % merge datasets
        [header,data,message_string] = RLW_merge_epochs(datasets,merge_idx);
        disp([message_string{1} ' of ' num2str(numel(datasets2merge)) ' datasets in ' seg_category{b} ' category.'])
        
        % create new name
        name_new = [seg_category{b} ' crop_120 ' base_name(1:end-4)];
        
        % save the data
        save([name_new '.mat'], 'data')
        
        % save the header
        header.name = name_new;
        save([name_new '.lw6'], 'header')
    end    
end
disp('>>>>>>> Done. <<<<<<<')

