%% GABA-AD: TEPs - TIME-FREQUENCY ANALYSIS
% Written by Dominika for GABA-AD project (2022)
% 
% Colection of scripts to perform TF decomposition of TEPs, visualize the
% outcome and run randomisation-based statistics on the results
% 
% Custom functions are included in the script, EEGLAB and letswave
% functions are being called from their directories
% 
% Output:
%   --> figures are saved in a folder 'GABA_YC_figures'
%   --> extracted values are saved in 'GABA_YC_variables'
% 
% 1) prepare single-trial TEP data
%       - load data and save them to a matlab structure (big!)
%       - downsample and save for letswave
%       - crop and export to the EEGLAB .set format
% 
% 2) short-time time-frequency decomposition
%       - 
% 
% ----- MEDICATION -----
% 3) subtrat pre-medciation dataset + pool ROIs 
%       -
% 
% 4) plot average change maps
%       -
% 
% 5) piared t-test 
%       -

% -------- SICI --------
% ) subtrat TS dataset + pool ROIs 
%       -
% 
% ) plot average change maps
%       -
% 
% ) piared t-test 
%       -

%% parameters
clear all; clc

% ----------- dataset -----------
study = 'GABA';
group = 'YC';
target = 'M1';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};
stimulus = {'CS' 'TS' 'ppTMS'};
n_chan = 30;
electrodes = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2',...
    'F7','F8','T7','T8','P7','P8','Fz','Cz','Pz','Iz','FC1','FC2',...
    'CP1','CP2','FC5','FC6','CP5','CP6','M1','M2','C1','C2'};
% --------------------------------

% get the starting filename prefix
prefix = sprintf('%s %s %s single_trial', study, group, target);

% navigate to the Git folder --> source of saved default files
folder_git = uigetdir(pwd, 'Choose the Git folder');

% visualization 
figure_counter = 1;

% input/output folders
folder_input = uigetdir(pwd, 'Choose the input folder');
folder_results = uigetdir(pwd, 'Choose the Results folder');
folder_figures = [folder_results '\GABA_' group '_figures'];
folder_output = [folder_results '\GABA_' group '_variables'];

%% 1) prepare single-trial TEP data
% ----- section input -----
window_extract = [-0.5, 0.5];
suffix = '500Hz';
% -------------------------
% create data output structure
if exist([folder_output '\GABA_YC_data.mat'])
    load([folder_output '\GABA_YC_data.mat'])
else
    GABA_YC_data = struct;
end

% loop through datasets
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            for s = 1:length(stimulus)
                % identify dataset
                if participant(p) < 10
                    subj = ['0' num2str(participant(p))];
                else
                    subj = num2str(participant(p));
                end
                filename = sprintf('%s %s %s %s %s', prefix, subj, medication{m}, time{t}, stimulus{s});

                % load data and lw header
                option = struct('filename', [folder_input '\' filename '.lw6']);
                lwdata = FLW_load.get_lwdata(option);
                
                % crop raw data data & save to output structure
                x_start = (window_extract(1) - lwdata.header.xstart)/lwdata.header.xstep;
                x_end = (window_extract(2) - lwdata.header.xstart)/lwdata.header.xstep;
                data_save = lwdata.data(:, 1:n_chan, :, :, :, x_start:x_end);
                data_save = permute(single(data_save),[1,2,6,3,4,5]);
                statement = ['GABA_YC_data.' target '.single_trial(m, t, s, p, 1:size(data_save, 1), 1:size(data_save, 2), 1:size(data_save, 3)) = data_save;'];
                eval(statement)
                save([folder_output '\GABA_YC_data.mat'], 'GABA_YC_data', '-v7.3')
                
                % downsample
                option = struct('x_dsratio', 4, 'suffix', suffix, 'is_save', 1);
                lwdata = FLW_downsample.get_lwdata(lwdata, option);
               
                % crop for EEGLAB
                x_start = (window_extract(1) - lwdata.header.xstart)/lwdata.header.xstep;
                x_end = (window_extract(2) - lwdata.header.xstart)/lwdata.header.xstep;
                lwdata.data = lwdata.data(:, 1:n_chan, :, :, :, x_start:x_end);
                data_extract = permute(single(lwdata.data),[2,6,1,3,4,5]);
                
                % modify header
                lwdata.header.xstart = window_extract(1);

                % create an EEGLAB structure & export
                load([folder_git '\EEGLAB_example.mat']);
                EEG = lw2eeglab(EEG, data_extract, lwdata.header);
                save([folder_results '\' study '_' group '_export\EEGLAB\' filename '.set'], 'EEG')
            end
        end
    end
    % update 
    fprintf('Subject n.%d finished.', participant(p))
end

% create subject average data & save
GABA_YC_data.M1.individual(:, :, :, :, :, :) = squeeze(mean(GABA_YC_data.M1.single_trial, 5));
save([folder_output '\GABA_YC_data.mat'], 'GABA_YC_data', '-v7.3')
clear m t s p subj filename lwdata data_save data_extract EEG x_start x_end statement

% update prefix
prefix = [suffix ' ' prefix];
clear suffix window_extract

%% 2) short-time time-frequency decomposition
% ----- section input -----
hanning = 0.25;
slide = 1;
freqs = [1, 45, 1];
output = 'power';
window_crop = [-0.5, 0.5];
baseline = [-0.5, -0.1];
suffix = {'stft' 'crop' 'full_z' 'avg' 'bl_z'};
% -------------------------
% loop through datasets
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            for s = 1:length(stimulus)
                % identify dataset
                if participant(p) < 10
                    subj = ['0' num2str(participant(p))];
                else
                    subj = num2str(participant(p));
                end
                filename = sprintf('%s %s %s %s %s', prefix, subj, medication{m}, time{t}, stimulus{s});
                disp(sprintf('Subject n.%d: %s, %s medication, %s stimulus', participant(p), medication{m}, time{t}, stimulus{s}))
                
                % load the dataset
                option = struct('filename', [folder_input '\' filename '.lw6']);                
                lwdata = FLW_load.get_lwdata(option);                
                
                % remove the target channel if it exists
                option = struct('type','channel', 'items', {electrodes}, 'suffix', '', 'is_save', 0);
                lwdata = FLW_selection.get_lwdata(lwdata, option);
                
                % STFT using FFT and Hanning window taper
                option = struct('hanning_width', hanning, 'sliding_step', slide, ...
                    'low_frequency', freqs(1), 'high_frequency', freqs(2), ...
                    'num_frequency_lines', (freqs(2) - freqs(1) + 1)/freqs(3), ...
                    'output', output, 'show_progress', 0, 'suffix', suffix{1}, 'is_save', 0);
                lwdata = FLW_STFT.get_lwdata(lwdata, option);
                fprintf('...STFT')
                
                % crop data
                option = struct('xcrop_chk', 1, 'xstart', window_crop(1), 'xend', window_crop(2), ...
                    'suffix', suffix{2}, 'is_save', 0);
                lwdata = FLW_crop_epochs.get_lwdata(lwdata, option);
                
                % compute z-scores based on the full-length epoch
                option = struct('operation', 'zscore', 'xstart', window_crop(1), 'xend', window_crop(2), ...
                    'suffix', suffix{3}, 'is_save', 0);
                lwdata = FLW_baseline.get_lwdata(lwdata, option);
                fprintf('...full-length correction')
                
                % average across trials
                option = struct('operation', 'average', 'suffix', suffix{4}, 'is_save', 0);
                lwdata = FLW_average_epochs.get_lwdata(lwdata, option);
                
                % baseline correct
                option = struct('operation', 'zscore', 'xstart', baseline(1), 'xend', baseline(2),...
                    'suffix', suffix{5}, 'is_save', 1);
                lwdata = FLW_baseline.get_lwdata(lwdata, option);
                fprintf('...baseline correction')
            end
        end
    end
end
clear p m t s subj filename option lwdata hanning slide freqs output window_crop baseline

% update prefix
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear s suffix 

%% 3) MEDICATION: subtrat pre-medciation dataset + pool ROIs 
% ----- section input -----
ROI(1).area = 'frontal'; ROI(2).area = 'central'; ROI(3).area = 'left_temporal'; ROI(4).area = 'right_temporal'; ROI(5).area = 'occipital'; 
ROI(1).electrodes = {'Fp1' 'Fp2' 'Fz' 'F3' 'F4' 'F7' 'F8'};
ROI(2).electrodes = {'FC1' 'FC2' 'Cz' 'C1' 'C2' 'CP1' 'CP2'};
ROI(3).electrodes = {'FC5' 'T7' 'C3' 'CP5'};
ROI(4).electrodes = {'FC6' 'T8' 'C4' 'CP6'};
ROI(5).electrodes = {'P3' 'P4' 'Pz' 'P7' 'P8' 'O1' 'O2' 'Iz'};
suffix = {'diff' 'merged_subj'};
% -------------------------
% subtract baseline from post-medication dataset & create ROIs
for s = 1:2
    for m = 1:length(medication) 
        disp(sprintf('stimulus %s: %s', stimulus{s}, medication{m}))
        for p = 1:length(participant)
            % identify datasets
            if participant(p) < 10
                subj = ['0' num2str(participant(p))];
            else
                subj = num2str(participant(p));
            end
            filename = sprintf('%s %s %s', prefix, subj, medication{m});

            % load the datasets 
            option = struct('filename',{{sprintf('%s\\%s %s %s.lw6', folder_input, filename, time{1}, stimulus{s}), ...
                sprintf('%s\\%s %s %s.lw6', folder_input, filename, time{2}, stimulus{s})}});
            set = FLW_load.get_lwdataset(option);
            
            % subtract baseline data from post-medication data
            option = struct('operation', 'sub', 'ref_dataset', 1, 'suffix', suffix{1}, 'is_save', 0);
            set = FLW_math.get_lwdataset(set,option);
            set = set(2);  
            
            % pool into ROIs and average channel
            for r = 1:length(ROI) + 1
                if r <= length(ROI)
                    % identify electrode locations
                    labels = {set.header.chanlocs.labels};
                    for e = 1:length(ROI(r).electrodes)
                        epos(e) = find(contains(labels, ROI(r).electrodes{e}) & cellfun('length', labels) == length(ROI(r).electrodes{e}));
                    end
                    
                    % pool selected electrodes
                    set.data(1, size(set.data, 2) + 1, 1, 1, :, :) = mean(set.data(1, epos, 1, 1, :, :), 2);
                    
                    % update header
                    set.header.chanlocs(length(set.header.chanlocs) + 1).labels = ROI(r).area;  
                else
                    % pool all electrodes
                    set.data(1, size(set.data, 2) + 1, 1, 1, :, :) = mean(set.data(1, :, 1, 1, :, :), 2);
                    
                    % update header
                    set.header.chanlocs(length(set.header.chanlocs) + 1).labels = 'average';  
                end
            end
            set.header.datasize(2) = size(set.data, 2);
            
            % append to overall dataset, change name
            lwdataset(p) = set;  
            lwdataset(p).header.name = sprintf('%s bl_z avg full_z crop stft 500Hz GABA YC M1 %s %s', suffix{1}, medication{m}, stimulus{s});
            
            % update
            fprintf('...%d', p)
        end
        
        % merge subjects to a single letswave dataset
        option = struct('type', 'epoch', 'suffix', suffix{2}, 'is_save', 1);
        lwdata = FLW_merge.get_lwdata(lwdataset, option);      
    end
end
clear s m p r e subj filename lwdataset lwdata set epos 

% update prefix
prefix = 'bl_z avg full_z crop stft 500Hz GABA YC M1'
for s = 1:length(suffix)
    prefix = [suffix{s} ' ' prefix];
end
clear s suffix

%% 4) MEDICATION: plot average change maps
% ----- section input -----
channel = 'average';
axlim = {[-200, 500], [1, 45]};      % limits of x (time) and y (frequency)
clim = {[-2 2], [-2 2]};          % limits of the colorscale {CS, TS}
% -------------------------
% loop through datasets
for s = 1:2
    for m = 1:length(medication)
        % load the data
        filename = sprintf('%s %s %s', prefix, medication{m}, stimulus{s});
        option = struct('filename', [folder_input '\' filename '.lw6']);                
        lwdata = FLW_load.get_lwdata(option);
                
        % plot the TF map
        fig = figure(figure_counter);
        plot_map(lwdata, channel, 'axlim', axlim, 'clim', clim{s});
        
        % name and save figure
        figure_name = sprintf('STFT_change_%s_%s', stimulus{s}, medication{m});
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')
        
        % update the counter
        figure_counter = figure_counter + 1; 
    end
end
clear s m channel axlim clim filename option lwdata fig figure_name

%% 5) MEDICATION: piared t-test
% ----- section input -----
alpha = 0.05;
alpha_cl = 0.05;
n_perm = 2000;
suffix = 'ttest';
% -------------------------
for s = 1:2
    % load datasets
    disp(sprintf('stimulus %s', stimulus{s}))
    option = struct('filename',{{sprintf('%s %s %s', prefix, medication{1}, stimulus{s}), ...
        sprintf('%s %s %s', prefix, medication{2}, stimulus{s})}});
    lwdataset = FLW_load.get_lwdataset(option);
    
    option = struct('test_type', 'paired sample', 'tails', 'both', 'ref_dataset', 1, 'alpha', alpha, ...
        'permutation', 1, 'cluster_threshold', alpha_cl,'num_permutations', n_perm, 'show_progress', 0, ...
        'multiple_sensor', 1, 'chan_dist', 0.27, 'suffix', suffix, 'is_save', 1);
    lwdataset = FLW_ttest.get_lwdataset(lwdataset, option);
end

% update prefix
prefix = [suffix ' ' prefix];
clear s option lwdataset alpha alpha_cl n_perm suffix

%% functions
function EEG = lw2eeglab(EEG, data, header)
    % deal with datasize
    EEG.setname = header.name;
    EEG.filename = [];
    EEG.filepath = [];
    EEG.nbchan = size(data, 1);
    EEG.trials = size(data, 3);
    EEG.pnts = size(data, 2);
    EEG.srate = 1/header.xstep;
    EEG.times = header.xstart+(0:EEG.pnts-1)*header.xstep;
    EEG.xmin = EEG.times(1);
    EEG.xmax = EEG.times(end);    
    EEG.data = data; 
    
    % deal with channel locations
    EEG.chanlocs = header.chanlocs(1:size(data, 1));
    EEG.chanlocs = rmfield(EEG.chanlocs, {'SEEG_enabled' 'topo_enabled', 'sph_theta_besa', 'sph_phi_besa'});
    EEG.chanlocs = orderfields(EEG.chanlocs, [1, 2, 3, 6, 7, 8, 4, 5]);
    [EEG.chanlocs.sph_radius] = deal(85);
    [EEG.chanlocs.type] = deal('EEG');
    C = num2cell(1:length(EEG.chanlocs));
    [EEG.chanlocs.urchan] = C{:};
    [EEG.chanlocs.ref] = deal([]);
    EEG.urchanlocs = EEG.chanlocs;
    
    % deal with events
    EEG.event = header.events;
    if ~isempty(EEG.event)
        [EEG.event.type] = EEG.event.code;
        EEG.event = rmfield(EEG.event,'code');
        temp = num2cell([EEG.event.latency]/header.xstep);
        [EEG.event.latency] = deal(temp{:});     
        [EEG.event.urevent] = EEG.event.epoch;
        EEG.event = orderfields(EEG.event, [3, 1, 4, 2]);
        EEG.urevent = EEG.event;
    end
end
function plot_map(lwdata, channel, varargin)
    % prepare colormap 
    map = colorcet('D1','N', 201);
    
    % prepare axes
    x = lwdata.header.xstart*1000 : lwdata.header.xstep*1000 : ((lwdata.header.datasize(6)-1)*lwdata.header.xstep + lwdata.header.xstart)*1000;
    y = lwdata.header.ystart : lwdata.header.ystep : (lwdata.header.datasize(5)-1)*lwdata.header.ystep + lwdata.header.ystart;
    
    % prepare data
    if strcmp(channel, 'average')
        data = squeeze(mean(lwdata.data, [1, 2]));
    else
        data = squeeze(mean(lwdata.data(:, channel, :, :, :, :), 1));
    end
    
    % launch the figure
    set(gcf, 'units','centimeters','position',[10 10 15 10], 'color', 'w');
    hold on
    
    % plot the map
    imagesc(x,y,data)
    
    % set colorscale
    colorbar
    colormap(map)
    c = find(strcmpi(varargin, 'clim'));
    if ~isempty(c)
        clim = varargin{c + 1};
        set(gca, 'clim', clim)
    end 
    
    % set x and y limits
    a = find(strcmpi(varargin, 'axlim'));
    if ~isempty(a)
        axlim = varargin{a + 1};
        xlim(axlim{1})
        ylim(axlim{2})
    end
    
    % add TMS stimulus
    line([0, 0], axlim{2}, 'Color', [0 0 0], 'LineWidth', 2.5, 'LineStyle', '--')
    
    % other parameters
    xlabel('time (ms)')
    ylabel('frequency (Hz)')
    set(gca, 'FontSize', 16) 
    set(gca, 'Layer', 'Top')
end


