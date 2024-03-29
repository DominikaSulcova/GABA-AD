%% GABA-AD: EXPORT + GROUP ANALYSIS
% Written by Dominika for GABA-AD project (2021)
%
% Colection of scripts to perform analysis of overall data: 
%   --> figures are saved in a folder 'GABA_YC_figures'
%   --> output variables are saved in a folder 'GABA_YC_statistics'
% 
% ----- DATA EXTRACTION -----
% Extracts subject info, parameters and outcome measures from all available 
% datasets and saves it in structure 'GABA_YC_results'
% 1) subject info
% 2) RS-EEG
% 3) M1 TEPs
% 4) MEPs
% 5) SICI
% 6) AG TEPs
% 7) export data for R
% 
% ----- RANOVA -----
% 8) arousal
% 9) rMT
% 
% ----- CORRELATIONS -----
% 14) TEP x MEP
% 15) TEP change x RS-EEG change
% 16) TEP SICI x MEP SICI 
% 17) TEP-peak SICI x MEP SICI 

%% parameters
clear all, clc

% ----- adjustable parameters -----
% dataset
group = 'YC';
participant = 1:20;
peaks_M1 = {'N17' 'P30' 'N45' 'P60' 'N100' 'P180'};
peaks_AG = {'P25' 'N40' 'P50' 'N75' 'N100' 'P180'};

% visualization 
time_window = [-0.05, 0.3];
shade = 0.2;

% statistics
alpha = 0.05;
z = 1.96;
% --------------------------------

% navigate to the Git folder --> source of saved default files
folder_git = uigetdir(pwd, 'Choose the Git folder');

% load default header
load([folder_git '\GABA_header_default.mat'])
for l = 1:size(header.chanlocs, 2)
    labels{l} = header.chanlocs(l).labels;
end
clear l

% visualization calculated params
figure_counter = 1;
xstep = header.xstep; 
xstart = header.xstart;
x = [time_window(1):xstep:time_window(2)];

% check for colour scheme
answer = questdlg('Do you want to choose a new colour scheme?', 'Colour scheme', 'YES', 'NO', 'NO'); 
switch answer
    case 'YES'
        a = 1;
        for p = 1:6
            for c = 1:length(current)
               colours(a, :) = uisetcolor; 
               a = a + 1;
            end
        end
        save('colours.mat', 'colours'); 
    case 'NO'
    load([folder_git '\GABA_colours.mat'])
end
clear a answer

% create output folders
folder_results = uigetdir(pwd, 'Choose the Results folder');
folder_input = [folder_results '\GABA_' group '_variables'];
folder_output = [folder_results '\GABA_' group '_statistics'];
if ~exist(folder_output)
    mkdir(folder_output)
end
folder_figures = [folder_results '\GABA_' group '_figures'];
if ~exist(folder_figures)
    mkdir(folder_figures)  
end

% TEPs
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};
stimulus = {'CS' 'TS' 'ppTMS'};
peaks = {'N17' 'P30' 'N45' 'P60' 'N100' 'P180'};

% RS-EEG
condition = {'open' 'closed'};
fband_names = {'delta' 'theta' 'alpha1' 'alpha2' 'alpha3' 'beta1' 'beta2' 'gamma'};
ROI = struct;
ROI(1).area = 'frontal'; ROI(2).area = 'central'; ROI(3).area = 'left_temporal'; ROI(4).area = 'right_temporal'; ROI(5).area = 'occipital'; 
ROI(1).electrodes = {'Fp1' 'Fp2' 'Fz' 'F3' 'F4' 'F7' 'F8'};
ROI(2).electrodes = {'FC1' 'FC2' 'Cz' 'C1' 'C2' 'CP1' 'CP2'};
ROI(3).electrodes = {'FC5' 'T7' 'C3' 'CP5'};
ROI(4).electrodes = {'FC6' 'T8' 'C4' 'CP6'};
ROI(5).electrodes = {'P3' 'P4' 'Pz' 'P7' 'P8' 'O1' 'O2' 'Iz'};

%% 1) prepare outcome structure, subject info
% load individual params
params = readtable([folder_input '\YC_parameters.csv']);

% fill in basic info
if ~exist([folder_output '\GABA_YC_results.mat'])
    % insert data from parameters
    GABA_YC_results = struct;
    for p = 1:length(participant)
        % subject info
        GABA_YC_results.subjects.age(p) = params{p, 'Age'}; 
        GABA_YC_results.subjects.sex(p) = params{p, 'Sex'};                         % 1 = female; 0 = male
        GABA_YC_results.subjects.handedness(p) = params{p, 'Handedness'};           % 1 = right; 0 = left
        GABA_YC_results.subjects.APOE(p) = params{p, 'APOE'};                       % 1 = homo/hetero eta4; 0 = no eta4

        % resting motor threshold
        GABA_YC_results.rmt(1).pre(p) = params{p, 'zyrtec_rMT_pre'}; 
        GABA_YC_results.rmt(1).post(p) = params{p, 'zyrtec_rMT_post'}; 
        GABA_YC_results.rmt(1).change(p) = params{p, 'zyrtec_rMT_post'} - params{p, 'xanax_rMT_pre'}; 
        GABA_YC_results.rmt(2).pre(p) = params{p, 'xanax_rMT_pre'}; 
        GABA_YC_results.rmt(2).post(p) = params{p, 'xanax_rMT_post'}; 
        GABA_YC_results.rmt(2).change(p) = params{p, 'xanax_rMT_post'} - params{p, 'xanax_rMT_pre'}; 
        
        % arousal        
        GABA_YC_results.arousal(1).assessment_1(p) = params{p, 14 + 1};
        GABA_YC_results.arousal(1).assessment_2(p) = params{p, 14 + 2};
        GABA_YC_results.arousal(1).assessment_3(p) = params{p, 14 + 3};
        GABA_YC_results.arousal(2).assessment_1(p) = params{p, 7 + 1};
        GABA_YC_results.arousal(2).assessment_2(p) = params{p, 7 + 2};
        GABA_YC_results.arousal(2).assessment_3(p) = params{p, 7 + 3};
    end
    
    % save the output variable
    save([folder_output '\GABA_YC_results.mat'], 'GABA_YC_results')
else
    load([folder_output '\GABA_YC_results.mat'])
end

%% 2) RS-EEG 
% load data
load([folder_input '\IAF.mat']) 
load([folder_input '\TF.mat']) 
load([folder_input '\fband.mat']) 
load([folder_input '\rsEEG_fband_amplitude.mat'])
load([folder_input '\AAC.mat']) 
load([folder_input '\AAC_change.mat'])
load([folder_input '\rsEEG_spect_exp.mat'])
load([folder_input '\beta.mat'])
load([folder_input '\beta_change.mat'])
load([folder_input '\delta.mat'])
load([folder_input '\delta_change.mat'])

% dataset info 
GABA_YC_results.rsEEG.ROI = ROI;
GABA_YC_results.rsEEG.fbands = fband_names; 
clear ROI

% IAF, TF, fband amplitude 
for p = 1:length(participant)
    for m = 1:length(medication)        
        % fband limits
        rows = (fband.subject == p & ...
            categorical(fband.medication) == medication{m});
        lim_i = fband{rows, 3:10};
        GABA_YC_results.rsEEG(m).limits(p, :) = lim_i;   
        
        % TF
        rows = (TF.subject == p & ...
            categorical(TF.medication) == medication{m});
        tf_i = TF.occipital(rows);
        GABA_YC_results.rsEEG(m).TF(p) = tf_i;        
        
        for t = 1:length(time)
            % IAF
            rows = (IAF.subject == p & ...
                categorical(IAF.medication) == medication{m} & ...
                categorical(IAF.time) == time{t});
            iaf_i = IAF{rows, 4:8};
            statement = ['GABA_YC_results.rsEEG(m).IAF.' time{t} '(p, :) = iaf_i;'];
            eval(statement)      
            
            % fband amplitude
            for c = 1:length(condition)
                for f = 1:length(fband_names)
                    rows = (fband_amplitude.subject == p & ...
                        categorical(fband_amplitude.medication) == medication{m} & ...
                        categorical(fband_amplitude.time) == time{t} & ...
                        categorical(fband_amplitude.condition) == condition{c} & ...
                        categorical(fband_amplitude.fband) == fband_names{f});
                    amp_i = fband_amplitude{rows, 6:10}';
                    statement = ['GABA_YC_results.rsEEG(m).amplitude.' time{t} '.' condition{c} '(p, :, f) = amp_i;'];
                    eval(statement)                
                end
            end 
        end
        
        % amplitude change
        for c = 1:length(condition)
            for f = 1:length(fband_names)                
                for t = 1:length(time)
                    statement = ['data_' time{t} ' = squeeze(GABA_YC_results.rsEEG(m).amplitude.' time{t} '.' condition{c} '(p, :, f));'];
                    eval(statement)
                end
                change_i = data_post - data_pre;

                % append
                statement = ['GABA_YC_results.rsEEG(m).amplitude.change.' condition{c} '(p, :, f) = change_i;'];
                eval(statement)
            end
        end         
    end
end
clear IAF iaf_i TF tf_i fband lim_i fband_amplitude amp_i change_i data_pre data_post rows statement f c 

% low beta - sigma peak
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            statement = ['GABA_YC_results.rsEEG(m).sigma.' time{t} '(p) = squeeze(beta(m, t, p, 1));'];
            eval(statement)            
        end
        
        % change = normalized
        GABA_YC_results.rsEEG(m).sigma.change(p) = squeeze(beta_change(m, p, 1));     
    end
end
clear beta beta_change p m t

% delta
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            statement = ['GABA_YC_results.rsEEG(m).delta.' time{t} '(p) = squeeze(delta(m, t, p));'];
            eval(statement)            
        end
        
        % change = normalized
        GABA_YC_results.rsEEG(m).delta.change(p) = squeeze(delta_change(m, p));     
    end
end
clear delta delta_change p m t

% alpha attenuation coeficient - broadband
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            % AAC
            statement = ['GABA_YC_results.rsEEG(m).AAC.' time{t} '(p) = squeeze(AAC(m, t, 4, p));'];
            eval(statement)            
        end
        
        % AAC change
        GABA_YC_results.rsEEG(m).AAC.change(p) = squeeze(AAC_change(m, 4, p));     
    end
end
clear AAC AAC_change statement p m t

% spectral exponent  
for p = 1:length(participant)
    for m = 1:length(medication)
        GABA_YC_results.rsEEG(m).SE.bands = {'low' 'high' 'broad'};        
        for t = 1:length(time)
            % SE
            for c = 1:length(condition)
                for a = 1:3
                    SE_i(a) = spect_exp(a).result.slope(m, t, c, p);
                end
                
                % append
                statement = ['GABA_YC_results.rsEEG(m).SE.' time{t} '.' condition{c} '(p, :) = SE_i;'];
                eval(statement)      
            end
        end
        
        % SE change
        for c = 1:length(condition)            
            for a = 1:3
                SE_change_i(a) = spect_exp(a).result.slopechange(m, c, p);
            end

            % append
            statement = ['GABA_YC_results.rsEEG(m).SE.change.' condition{c} '(p, :) = SE_change_i;'];
            eval(statement)      
        end    
        
        % SE change in %
        for c = 1:length(condition)   
            for a = 1:3
                % extract data
                statement = ['SE_pre = GABA_YC_results.rsEEG(m).SE.pre.' condition{c} '(p, a);'];
                eval(statement)
                statement = ['SE_post = GABA_YC_results.rsEEG(m).SE.post.' condition{c} '(p, a);'];
                eval(statement)

                % calculate change and append
                statement = ['GABA_YC_results.rsEEG(m).SE.change_percent.' condition{c} '(p, a) = SE_post/SE_pre*100 - 100;'];
                eval(statement)     
            end
        end   
    end
end
clear a m t c spect_exp SE_i SE_change_i data_pre statement SE_pre SE_post

% save
save([folder_output '\GABA_YC_results.mat'], 'GABA_YC_results', '-append')

%% 3) M1 TEPs
% load data
load([folder_input '\epochs.mat']);
load([folder_input '\GABA_YC_M1_TEPs.mat'], 'GABA_TEP_peaks');

% fill in TEP related info
for p = 1:length(participant)
    % placebo
    GABA_YC_results.TEP_M1(1).ICA.pre(p) = params{p, 'zyrtec_ICA_pre'};
    GABA_YC_results.TEP_M1(1).ICA.post(p) = params{p, 'zyrtec_ICA_post'};
    GABA_YC_results.TEP_M1(1).epochs.pre(p, :) = squeeze(epochs(p, 1, 1, :))';
    GABA_YC_results.TEP_M1(1).epochs.post(p, :) = squeeze(epochs(p, 1, 2, :))';
    
    % alprazolam
    GABA_YC_results.TEP_M1(2).ICA.pre(p) = params{p, 'xanax_ICA_pre'};
    GABA_YC_results.TEP_M1(2).ICA.post(p) = params{p, 'xanax_ICA_post'};
    GABA_YC_results.TEP_M1(2).epochs.pre(p, :) = squeeze(epochs(p, 2, 1, :))';
    GABA_YC_results.TEP_M1(2).epochs.post(p, :) = squeeze(epochs(p, 2, 2, :))';
end
clear params epochs

% fill in peak amplitude and latency
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            % peak amplitude
            statement = ['GABA_YC_results.TEP_M1(m).amplitude.' time{t} '(p, :, :) = squeeze(GABA_TEP_peaks.amplitude_peak(m, t, :, p, :));'];
            eval(statement)
            
            % peak amplitude
            statement = ['GABA_YC_results.TEP_M1(m).latency.' time{t} '(p, :, :) = squeeze(GABA_TEP_peaks.latency(m, t, :, p, :));'];
            eval(statement)
        end
        
        % calculate change 
        for s = 1:length(stimulus)
            for k = 1:length(peaks_M1)
                % define polarity
                if peaks_M1{k}(1) == 'P'
                    polarity = 1;
                elseif peaks_M1{k}(1) == 'N'
                    polarity = -1;
                end

                % change in amplitude in ?V
                amp_pre = GABA_YC_results.TEP_M1(m).amplitude.pre(p, s, k);
                amp_post = GABA_YC_results.TEP_M1(m).amplitude.post(p, s, k);
                GABA_YC_results.TEP_M1(m).amplitude.change(p, s, k) = (amp_post - amp_pre) * polarity;
                               
                % change in latency in ms
                lat_pre = GABA_YC_results.TEP_M1(m).latency.pre(p, s, k);
                lat_post = GABA_YC_results.TEP_M1(m).latency.post(p, s, k);
                GABA_YC_results.TEP_M1(m).latency.change(p, s, k) = (lat_post - lat_pre);
            end   
        end
    end
end
clear p m t s k statement polarity GABA_TEP_peaks amp_pre amp_post lat_pre lat_post

% save
save([folder_output '\GABA_YC_results.mat'], 'GABA_YC_results', '-append')

%% 4) MEPs
% load data
load([folder_input '\GABA_MEP.mat']);

% fill in MEP related info
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            rows = (output.subject == p & ...
                categorical(output.medication) == medication{m} & ...
                categorical(output.time) == time{t});     
            
            % zero trials included
            statement = ['GABA_YC_results.MEP(m).epochs.' time{t} '(p, :) = output.epochs_zero(rows);'];
            eval(statement)
            
            % zero trials discarded
            statement = ['GABA_YC_results.MEP(m).epochs_nozero.' time{t} '(p, :) = output.epochs_nozero(rows);'];
            eval(statement)
        end
    end
end
clear p m t rows statement

% fill in peak amplitude and number of included trials --> zero response included
for p = 1:length(participant)
    for m = 1:length(medication)   
    % loop through the output table
        for t = 1:length(time)
            % choose the row
            rows = (output.subject == p & ...
                categorical(output.medication) == medication{m} & ...
                categorical(output.time) == time{t});

            % zero trials included
            statement = ['GABA_YC_results.MEP(m).amplitude.' time{t} '(p, :) = output.amplitude_zero(rows);'];
            eval(statement)
            
            % zero trials discarded
            statement = ['GABA_YC_results.MEP(m).amplitude_nozero.' time{t} '(p, :) = output.amplitude_nozero(rows);'];
            eval(statement)
                            
        end

        % MEP change - zero trials included
        rows_pre = (output.subject == p & ...
                categorical(output.medication) == medication{m} & ...
                categorical(output.time) == 'pre');
        rows_post = (output.subject == p & ...
                categorical(output.medication) == medication{m} & ...
                categorical(output.time) == 'post');
        GABA_YC_results.MEP(m).amplitude.change(p, :) = round(output.amplitude_zero(rows_post) ./ output.amplitude_zero(rows_pre) * 100, 2);     
        
        % MEP change - zero trials discarded
        rows_pre = (output.subject == p & ...
                categorical(output.medication) == medication{m} & ...
                categorical(output.time) == 'pre');
        rows_post = (output.subject == p & ...
                categorical(output.medication) == medication{m} & ...
                categorical(output.time) == 'post');
        GABA_YC_results.MEP(m).amplitude_nozero.change(p, :) = round(output.amplitude_nozero(rows_post) ./ output.amplitude_nozero(rows_pre) * 100, 2);   
    end
end
clear p m t output rows statement rows_pre rows_post 

% save
save([folder_output '\GABA_YC_results.mat'], 'GABA_YC_results', '-append')

%% 5) SICI 
% load data
load([folder_input '\GABA_YC_M1_TEPs.mat'], 'GABA_SICI_peaks', 'GABA_SICI');

% calculate change in MEPs
for p = 1:length(participant)
    for m = 1:length(medication) 
        % baseline
        GABA_YC_results.SICI(m).MEP.pre(p) = GABA_YC_results.MEP(m).amplitude.pre(p, 2)/GABA_YC_results.MEP(m).amplitude.pre(p, 1)*100;
        
        % post medication
        GABA_YC_results.SICI(m).MEP.post(p) = GABA_YC_results.MEP(m).amplitude.post(p, 2)/GABA_YC_results.MEP(m).amplitude.post(p, 1)*100;
    end
end
clear p m 

% calculate change in TEPs for each TEP peak
for p = 1:length(participant)
    for m = 1:length(medication) 
        for k = 1:6
            % define polarity
            if mod(k, 2) == 1
                polarity = -1;
            else
                polarity = 1;
            end
            
            % baseline
            GABA_YC_results.SICI(m).TEP.pre(p, k) = (GABA_YC_results.TEP_M1(m).amplitude.pre(p, 3, k) ...
                - GABA_YC_results.TEP_M1(m).amplitude.pre(p, 2, k)) * polarity;

            % post medication
            GABA_YC_results.SICI(m).TEP.post(p, k) = (GABA_YC_results.TEP_M1(m).amplitude.post(p, 3, k) ...
                - GABA_YC_results.TEP_M1(m).amplitude.post(p, 2, k)) * polarity;
        end
    end
end
clear p m k polarity

% fill in TEP SICI peak amplitude
for p = 1:length(participant)
    for m = 1:length(medication) 
        for k = [1 2]
        % baseline
        GABA_YC_results.SICI(m).TEP_SICI.pre(p, k) = GABA_SICI_peaks.amplitude_peak(m, 1, p, k);

        % post medication
        GABA_YC_results.SICI(m).TEP_SICI.post(p, k) = GABA_SICI_peaks.amplitude_peak(m, 2, p, k);
        
        % change
        if mod(k, 2) == 1
            polarity = -1;
        else
            polarity = 1;
        end
        GABA_YC_results.SICI(m).TEP_SICI.change(p, k) = (GABA_SICI_peaks.amplitude_peak(m, 2, p, k) ...
            - GABA_SICI_peaks.amplitude_peak(m, 1, p, k)) * polarity;
        end
    end
end
clear p m k polarity

% fill in TEP SICI GFP_AUC
for p = 1:length(participant)
    for m = 1:length(medication) 
        % baseline
        GABA_YC_results.SICI(m).GFP_AUC.pre(p) = GABA_SICI.GFP_AUC(m, 1, p);

        % post medication
        GABA_YC_results.SICI(m).GFP_AUC.post(p) = GABA_SICI.GFP_AUC(m, 2, p);

        % change
        GABA_YC_results.SICI(m).GFP_AUC.change(p) = (GABA_SICI.GFP_AUC(m, 2, p)/GABA_SICI.GFP_AUC(m, 1, p)) * 100;
    end
end
clear p m 

% save
save([folder_output '\GABA_YC_results.mat'], 'GABA_YC_results', '-append')
clear GABA_SICI_peaks GABA_SICI

%% 6) AG TEPs
% load data
load([folder_input '\GABA_YC_AG_TEPs.mat'], 'GABA_TEP_peaks');

% fill in peak amplitude and latency
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            % peak amplitude
            statement = ['GABA_YC_results.TEP_AG(m).amplitude.' time{t} '(p, :) = squeeze(GABA_TEP_peaks.amplitude_peak(m, t, p, :));'];
            eval(statement)
            
            % peak amplitude
            statement = ['GABA_YC_results.TEP_AG(m).latency.' time{t} '(p, :) = squeeze(GABA_TEP_peaks.latency(m, t, p, :));'];
            eval(statement)
        end
        
         % calculate change 
        for k = 1:length(peaks_AG)
            % define polarity
            if peaks_AG{k}(1) == 'P'
                polarity = 1;
            elseif peaks_AG{k}(1) == 'N'
                polarity = -1;
            end

            % change in amplitude in ?V
            amp_pre = GABA_YC_results.TEP_AG(m).amplitude.pre(p, k);
            amp_post = GABA_YC_results.TEP_AG(m).amplitude.post(p, k);
            GABA_YC_results.TEP_AG(m).amplitude.change(p, k) = (amp_post - amp_pre) * polarity;

            % change in latency in ms
            lat_pre = GABA_YC_results.TEP_AG(m).latency.pre(p, k);
            lat_post = GABA_YC_results.TEP_AG(m).latency.post(p, k);
            GABA_YC_results.TEP_AG(m).latency.change(p, k) = (lat_post - lat_pre); 
        end    
    end
end
clear p m t k statement polarity GABA_TEP_peaks amp_pre amp_post lat_pre lat_post

% save
save([folder_output '\GABA_YC_results.mat'], 'GABA_YC_results', '-append')

%% 7) export data for R
% update 
time_new = [time, {'change'}];
stimulus_new = {'M1_CS' 'M1_TS' 'AG'};

% ----- pharmacological activation of GABAARs: TEPs -----
% prepare empty table
table_TEP = table;

% cycle through entries (= rows)
row_cnt = 1;
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time_new)
            for s = 1:length(stimulus_new)
                for k = 1:6
                    % fill corresponding line:
                    % condition info
                    table_TEP.subject(row_cnt) = participant(p);             
                    table_TEP.medication(row_cnt) = medication(m);
                    table_TEP.time(row_cnt) = time_new(t);
                    table_TEP.stimulus(row_cnt) = stimulus_new(s);
                    if table_TEP.stimulus{row_cnt}(1) == 'M'
                        table_TEP.peak(row_cnt) = peaks_M1(k);
                    elseif table_TEP.stimulus{row_cnt}(1) == 'A'
                        table_TEP.peak(row_cnt) = peaks_AG(k);
                    end

                    % covariate - rMT
                    statement = ['table_TEP.rMT(row_cnt) = GABA_YC_results.rmt(m).' time_new{t} '(p);'];
                    eval(statement)
                    
                    % outcome variables 
                    if s == 1 || s == 2
                            % amplitude
                            statement = ['table_TEP.amplitude(row_cnt) = GABA_YC_results.TEP_M1(m).amplitude.' time_new{t} '(p, s, k);'];
                            eval(statement)
                            
                            % latency
                            statement = ['table_TEP.latency(row_cnt) = GABA_YC_results.TEP_M1(m).latency.' time_new{t} '(p, s, k);'];
                            eval(statement)
                            
                    elseif s == 3
                            % amplitude
                            statement = ['table_TEP.amplitude(row_cnt) = GABA_YC_results.TEP_AG(m).amplitude.' time_new{t} '(p, k);'];
                            eval(statement)
                            
                            % latency
                            statement = ['table_TEP.latency(row_cnt) = GABA_YC_results.TEP_AG(m).latency.' time_new{t} '(p, k);'];
                            eval(statement)
                    end

                    % update row count
                    row_cnt = row_cnt + 1;
                end
            end
        end
    end
end
clear row_cnt p m t s k statement 

% % for MATLAB2016 and lower
% GABA_YC_medication_TEP = table;
% GABA_YC_medication_TEP.subject = table_TEP.subject';
% GABA_YC_medication_TEP.medication = table_TEP.medication';
% GABA_YC_medication_TEP.time = table_TEP.time';
% GABA_YC_medication_TEP.stimulus = table_TEP.stimulus';
% GABA_YC_medication_TEP.peak = table_TEP.peak';
% GABA_YC_medication_TEP.rMT = table_TEP.rMT';
% GABA_YC_medication_TEP.amplitude = table_TEP.amplitude';
% GABA_YC_medication_TEP.latency = table_TEP.latency';

% for newer MATLAB
GABA_YC_medication_TEP = table_TEP;

% save as CSV
writetable(GABA_YC_medication_TEP, [folder_output '\GABA_YC_medication_TEP.csv'], 'Delimiter', ',');

% ----- pharmacological activation of GABAARs: MEPs -----
% prepare empty table
table_MEP = table;

% cycle through entries (= rows)
row_cnt = 1;
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time_new)
            % fill corresponding line:
            % condition info
            table_MEP.subject(row_cnt) = participant(p);             
            table_MEP.medication(row_cnt) = medication(m);
            table_MEP.time(row_cnt) = time_new(t);

            % covariate - rMT
            statement = ['table_MEP.rMT(row_cnt) = GABA_YC_results.rmt(m).' time_new{t} '(p);'];
            eval(statement)

            % outcome variable - amplitude
            statement = ['table_MEP.amplitude(row_cnt) = GABA_YC_results.MEP(m).amplitude.' time_new{t} '(p, 1);'];
            eval(statement)

            % update row count
            row_cnt = row_cnt + 1;
        end
    end
end
clear row_cnt p m t s k statement 

% % for MATLAB2016 and lower
% GABA_YC_medication_MEP = table;
% GABA_YC_medication_MEP.subject = table_MEP.subject';
% GABA_YC_medication_MEP.medication = table_MEP.medication';
% GABA_YC_medication_MEP.time = table_MEP.time';
% GABA_YC_medication_MEP.rMT = table_MEP.rMT';
% GABA_YC_medication_MEP.amplitude = table_MEP.amplitude';

% for newer MATLAB
GABA_YC_medication_MEP = table_MEP;

% save as CSV
writetable(GABA_YC_medication_MEP, [folder_output '\GABA_YC_medication_MEP.csv'], 'Delimiter', ',');

% ----- pharmacological activation of GABAARs: RS-EEG -----
% prepare empty table
table_rsEEG = table;

% cycle through entries (= rows)
row_cnt = 1;
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time_new)
            % fill corresponding line:
            % condition info
            table_rsEEG.subject(row_cnt) = participant(p);             
            table_rsEEG.medication(row_cnt) = medication(m);
            table_rsEEG.time(row_cnt) = time_new(t);

            % outcome variables - delta
            statement = ['table_rsEEG.delta(row_cnt) = GABA_YC_results.rsEEG(m).delta.' time_new{t} '(p);'];
            eval(statement)
            
            % outcome variables - sigma
            statement = ['table_rsEEG.sigma(row_cnt) = GABA_YC_results.rsEEG(m).sigma.' time_new{t} '(p);'];
            eval(statement)

            % outcome variables - AAC
            statement = ['table_rsEEG.AAC(row_cnt) = GABA_YC_results.rsEEG(m).AAC.' time_new{t} '(p);'];
            eval(statement)
            
            % outcome variables - SE
            statement = ['table_rsEEG.SE(row_cnt) = GABA_YC_results.rsEEG(m).SE.' time_new{t} '.closed(p, 1);'];
            eval(statement)
            
            % update row count
            row_cnt = row_cnt + 1;
        end
    end
end
clear row_cnt p m t statement 

% % for MATLAB2016 and lower
% GABA_YC_medication_rsEEG = table;
% GABA_YC_medication_rsEEG.subject = table_rsEEG.subject';
% GABA_YC_medication_rsEEG.medication = table_rsEEG.medication';
% GABA_YC_medication_rsEEG.time = table_rsEEG.time';
% GABA_YC_medication_rsEEG.delta = table_rsEEG.delta';
% GABA_YC_medication_rsEEG.sigma = table_rsEEG.sigma';
% GABA_YC_medication_rsEEG.AAC = table_rsEEG.AAC';
% GABA_YC_medication_rsEEG.SE = table_rsEEG.SE';

% for newer MATLAB
GABA_YC_medication_rsEEG = table_rsEEG;

% save as CSV
writetable(GABA_YC_medication_rsEEG, [folder_output '\GABA_YC_medication_rsEEG.csv'], 'Delimiter', ',');

% ----- paired-pulse activation of GABAARs: baseline SICI -----
% prepare empty table
table_SICI_bl = table;

% cycle through entries (= rows)
row_cnt = 1;
for p = 1:length(participant)
    for m = 1:length(medication)
        for s = 2:3;
            for k = 1:6
                % fill corresponding line:
                % condition info
                table_SICI_bl.subject(row_cnt) = participant(p);             
                table_SICI_bl.medication(row_cnt) = medication(m);
                table_SICI_bl.stimulus(row_cnt) = stimulus(s);
                table_SICI_bl.peak(row_cnt) = peaks_M1(k);

                % outcome variables - raw TEP amplitude 
                table_SICI_bl.amplitude(row_cnt) = GABA_YC_results.TEP_M1(m).amplitude.pre(p, s, k);
                
                % outcome variables - raw MEP amplitude 
                table_SICI_bl.amplitude_MEP(row_cnt) = GABA_YC_results.MEP(m).amplitude.pre(p, s-1);
                
                % update row count
                row_cnt = row_cnt + 1;
            end
        end
    end
end
clear row_cnt p m s k  

% % for MATLAB2016 and lower
% GABA_YC_SICI_bl = table;
% GABA_YC_SICI_bl.subject = table_SICI_bl.subject';
% GABA_YC_SICI_bl.medication = table_SICI_bl.medication';
% GABA_YC_SICI_bl.stimulus = table_SICI_bl.stimulus';
% GABA_YC_SICI_bl.peak = table_SICI_bl.peak';
% GABA_YC_SICI_bl.amplitude = table_SICI_bl.amplitude';

% for newer MATLAB
GABA_YC_SICI_bl = table_SICI_bl;

% save as CSV
writetable(GABA_YC_SICI_bl, [folder_output '\GABA_YC_SICI_baseline.csv'], 'Delimiter', ',');

% ----- paired-pulse activation of GABAARs: SICI x medication -----
% prepare empty table
table_SICI_med = table;

% cycle through entries (= rows)
row_cnt = 1;
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time);
            for k = 1:6
                % fill corresponding line:
                % condition info
                table_SICI_med.subject(row_cnt) = participant(p);             
                table_SICI_med.medication(row_cnt) = medication(m);
                table_SICI_med.time(row_cnt) = time(t);
                table_SICI_med.peak(row_cnt) = peaks_M1(k);

                % outcome variables - TEP SICI 
                statement = ['table_SICI_med.SICI(row_cnt) = GABA_YC_results.SICI(m).TEP.' time{t} '(p, k);'];
                eval(statement)
                
                % outcome variables - TEP SICI 
                statement = ['table_SICI_med.SICI_MEP(row_cnt) = GABA_YC_results.SICI(m).MEP.' time{t} '(p);'];
                eval(statement)
                
                % covariate - rMT
                statement = ['table_SICI_med.rMT(row_cnt) = GABA_YC_results.rmt(m).' time{t} '(p);'];
                eval(statement)
                
                % update row count
                row_cnt = row_cnt + 1;
            end
        end
    end
end
clear row_cnt p m t k statement  

% % for MATLAB2016 and lower
% GABA_YC_SICI_med = table;
% GABA_YC_SICI_med.subject = table_SICI_med.subject';
% GABA_YC_SICI_med.medication = table_SICI_med.medication';
% GABA_YC_SICI_med.time = table_SICI_med.time';
% GABA_YC_SICI_med.SICI = table_SICI_med.SICI';
% GABA_YC_SICI_med.rMT = table_SICI_med.rMT';

% for newer MATLAB
GABA_YC_SICI_med = table_SICI_med;

% save as CSV
writetable(GABA_YC_SICI_med, [folder_output '\GABA_YC_SICI_medication.csv'], 'Delimiter', ',');

clear table_TEP table_MEP table_rsEEG table_SICI_bl table_SICI_med 

%% 8) RANOVA: arousal
timepoint = {'1.5h' '2h' '2.5h'};

% extract data
for m = 1:length(medication)
    for p = 1:length(participant)
        for t = 1:length(timepoint)
            statement = ['data_arousal(m, p, t) = GABA_YC_results.arousal(m).assessment_' num2str(t) '(p);'];
            eval(statement)
        end
    end
end
clear m p t statement

% extract average values
for m = 1:length(medication)
    for t = 1:length(timepoint)
        data_arousal_avg(m, t, 1) = mean(data_arousal(m, :, t), 'omitnan');
        data_arousal_avg(m, t, 2) = std(data_arousal(m, :, t), 'omitnan');      % STD
        n = sum(~isnan(data_arousal(m, :, t)));
        data_arousal_avg(m, t, 3) = data_arousal_avg(m, t, 2)/sqrt(n);              % SEM
        data_arousal_avg(m, t, 4) = data_arousal_avg(m, t, 3) * z;                  % 95% confidence interval      
    end
end
clear t m n
% pd = fitdist(data_med(1, :, 1)','Normal')
% CI = paramci(pd);

% launch the figure
fig = figure(figure_counter);
hold on

% plot data with errorbars
x = 1:length(timepoint);    
for m = 1:length(medication)  
    % plot
    err(m) = errorbar(x, data_arousal_avg(m, :, 1), data_arousal_avg(m, :, 4));

    % add parameters
    err(m).Color = colours((m - 1)*2 + 2, :);
    err(m).LineWidth = 2;
    err(m).Marker = 'o';
    err(m).MarkerFaceColor = colours((m - 1)*2 + 2, :);
    err(m).MarkerSize = 10;
end

% set other parameters
fig_name = 'Subjective arousal assessment';
title(fig_name, 'FontWeight', 'bold', 'FontSize', 16)
set(gca, 'Fontsize', 16)
legend(err, medication, 'Location', 'southeast', 'fontsize', 16, 'EdgeColor', 'none')
xlabel('timepoint')
set(gca, 'xtick', 1:length(timepoint), 'xticklabel', timepoint)
ylim([0 10])
xlim([0.75 length(timepoint)+0.25])
ax = gca; ax.YGrid = 'on';
hold off
clear x ax err fig fig_name m 

% save figure
fig_name = 'arousal';
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name  '.png'])

% update counter
figure_counter = figure_counter + 1;

% ----- test for differences - rm-ANOVA -----
% prepare data table
data_table = table; data_table.participant = participant'; 
for m = 1:length(medication)
    for t = 1:length(timepoint)
        statement = ['data_table.' medication{m} '_t' num2str(t) ' = data_arousal(m, :, t)'';'];
        eval(statement)
    end
end
clear m t statement 
 
% within-subjects design
wd = table([repmat(medication(1), length(timepoint), 1); repmat(medication(2), length(timepoint), 1)],...
    [1:length(timepoint) 1:length(timepoint)]','VariableNames',{'treatment','timepoint'});
wd.treatment = categorical(wd.treatment); wd.timepoint = categorical(wd.timepoint);

% rm-ANOVA
rm = fitrm(data_table, 'placebo_t1-alprazolam_t3 ~ 1', 'WithinDesign', wd);

% check for approx. normal distribution
figure(figure_counter)
histogram(reshape(data_table{:,2:end},[],1), 10)

% check for sphericity
rm.mauchly

% results table
ranova_arousal = ranova(rm, 'WithinModel', 'treatment*timepoint');

% plot as boxplot
fig = figure(figure_counter);
boxplot(data_table{:,2:end})
set(gca, 'FontSize', 14)
set(gca, 'xtick', 1:6, 'xticklabel', [1 2 3 1 2 3])
text(1.5, -0.5, 'placebo', 'FontSize', 14)
text(4.5, -0.5, 'alprazolam', 'FontSize', 14)
ylabel('arousal assessment')

% save figure
fig_name = 'arousal_boxplot';
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name  '.png'])

% update counter
figure_counter = figure_counter + 1;

% save ANOVA output
ranova_arousal = ranova_arousal(:, 1:6);
writetable(ranova_arousal, [folder_output '\GABA_YC_ranova.xlsx'],  'Sheet', 'arousal')

% add mean values
writematrix(data_arousal_avg(:, :, 1), [folder_output '\GABA_YC_ranova.xlsx'],...
    'Sheet', 'arousal', 'WriteMode', 'append')
writematrix(data_arousal_avg(:, :, 3), [folder_output '\GABA_YC_ranova.xlsx'],...
    'Sheet', 'arousal', 'WriteMode', 'append')
clear timepoint data_arousal data_arousal_avg data_table fig fig_name rm wd ranova_arousal

%% 9) RANOVA: rMT 
% ----- extract data -----
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            statement = ['data_rMT(p, m, t) = GABA_YC_results.rmt(m).' time{t} '(p);'];
            eval(statement)
        end
    end
end
clear m p statement

% extract average values
for m = 1:length(medication)
    for t = 1:length(time)
        data_rMT_avg(m, t, 1) = mean(data_rMT(:, m, t), 'omitnan');
        data_rMT_avg(m, t, 2) = std(data_rMT(:, m, t), 'omitnan');      % STD
        n = sum(~isnan(data_rMT(:, m, t)));
        data_rMT_avg(m, t, 3) = data_rMT_avg(m, t, 2)/sqrt(n);              % SEM
        data_rMT_avg(m, t, 4) = data_rMT_avg(m, t, 3) * z;                  % 95% confidence interval      
    end
end
clear t m n
% pd = fitdist(data_med(1, :, 1)','Normal')
% CI = paramci(pd);

% ----- plot boxplots with lines -----
% prepapare data
for m = 1:length(medication)
    for t = 1:length(time)
        data_visual(:, (m-1)*2 + t) = data_rMT(:, m, t);
    end
end
clear m t

% plot group boxplot
fig = figure(figure_counter);     
hold on
boxplot(data_visual, 'color', colours)

% plot the lines - placebo
for p = 1:length(participant)
    p_placebo(p) = plot([1 2], data_visual(p, [1 2]), '-o',...
        'Color', [0.75, 0.75, 0.75],...
        'MarkerSize', 10,...
        'MArkerEdge', 'none');
    hold on
end

% plot the lines - alprazolam
for p = 1:length(participant)
    p_alprazolam(p) = plot([3 4], data_visual(p, [3 4]), '-',...
        'Color', [0.75, 0.75, 0.75],...
        'MarkerSize', 10,...
        'MArkerEdge', 'none');
    hold on
end

% plot the markers
for b = 1:size(data_visual, 2)
    scat(b) = scatter(repelem(b, length(participant)), data_visual(:, b),...
        75, colours(b, :), 'filled');
    hold on
end

% add parameters
figure_title = 'RESTING MOTOR THRESHOLD';
yl = get(gca, 'ylim');
set(gca, 'xtick', 1:4, 'xticklabel', {'pre' 'post' 'pre' 'post'})
set(gca, 'Fontsize', 14)
title(figure_title, 'FontWeight', 'bold', 'FontSize', 14)
xlabel('time relative to medication'); ylabel('rMT (%MSO)');
ylim([yl(1), yl(2) + (yl(2) - yl(1))*0.15])
hold on

% add text
txt(1) = text(1.2, yl(2) + (yl(2) - yl(1))*0.05, 'placebo', 'fontsize', 14, 'color', colours(2, :));
txt(2) = text(3.2, yl(2) + (yl(2) - yl(1))*0.05, 'alprazolam', 'fontsize', 14, 'color', colours(4, :));
hold off

% save figure
fig_name = 'rMT_boxplot';
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name  '.png'])

% update the counter
figure_counter = figure_counter + 1;   
clear data_visual b fig fig_name figure_title p p_alprazolam p_placebo scat txt yl
    
% ----- test for differences - rm-ANOVA -----
% prepare data table
data_table = table;
for m = 1:length(medication)
    for t = 1:length(time)
        statement = ['data_table.' medication{m} '_' time{t} ' = data_rMT(:, m, t);'];
        eval(statement)
    end
end
clear m t statement 

% within-subjects design
wd = table([repmat(medication(1), length(time), 1); repmat(medication(2), length(time), 1)],...
    [1:length(time) 1:length(time)]','VariableNames',{'medication','time'});
wd.medication = categorical(wd.medication); wd.time = categorical(wd.time);

wd = table([repmat(medication(1), 1, 1); repmat(medication(2), 1, 1)],'VariableNames',{'medication'});
wd.medication = categorical(wd.medication);
rm = fitrm(data_table, 'placebo_pre-alprazolam_pre ~ 1', 'WithinDesign', wd);
ranova_rMT = ranova(rm, 'WithinModel', 'medication');

% rm-ANOVA
rm = fitrm(data_table, 'placebo_pre-alprazolam_post ~ 1', 'WithinDesign', wd);

% check for approx. normal distribution
figure(figure_counter)
histogram(reshape(data_table{:,2:end},[],1), 10)

% check for sphericity
rm.mauchly

% results table
ranova_rMT = ranova(rm, 'WithinModel', 'medication*time');

% save ANOVA output 
ranova_rMT = ranova_rMT(:, 1:6);
writetable(ranova_rMT, [folder_output '\GABA_YC_ranova.xlsx'],  'Sheet', 'rMT')

% add mean values
writematrix(data_rMT_avg(:, :, 1), [folder_output '\GABA_YC_ranova.xlsx'],...
    'Sheet', 'rMT', 'WriteMode', 'append')
writematrix(data_rMT_avg(:, :, 3), [folder_output '\GABA_YC_ranova.xlsx'],...
    'Sheet', 'rMT', 'WriteMode', 'append')
clear data_rMT data_rMT_avg data_table rm wd ranova_rMT

%% 10) RANOVA: rsEEG - sigma
% ----- extract data -----
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            statement = ['data(p, m, t) = GABA_YC_results.rsEEG(m).sigma.' time{t} '(p);'];
            eval(statement)
        end
    end
end
clear m p t statement

% extract average values
for m = 1:length(medication)
    for t = 1:length(time)
        data_avg(m, t, 1) = mean(data(:, m, t), 'omitnan');
        data_avg(m, t, 2) = std(data(:, m, t), 'omitnan');      % STD
        n = sum(~isnan(data(:, m, t)));
        data_avg(m, t, 3) = data_avg(m, t, 2)/sqrt(n);              % SEM
        data_avg(m, t, 4) = data_avg(m, t, 3) * z;                  % 95% confidence interval      
    end
end
clear t m n
% pd = fitdist(data_med(1, :, 1)','Normal')
% CI = paramci(pd);
    
% ----- test for differences - rm-ANOVA -----
% prepare data table
data_table = table;
for m = 1:length(medication)
    for t = 1:length(time)
        statement = ['data_table.' medication{m} '_' time{t} ' = data(:, m, t);'];
        eval(statement)
    end
end
clear m t statement 

% within-subjects design
wd = table([repmat(medication(1), length(time), 1); repmat(medication(2), length(time), 1)],...
    [1:length(time) 1:length(time)]','VariableNames',{'medication','time'});
wd.medication = categorical(wd.medication); wd.time = categorical(wd.time);

% rm-ANOVA
rm = fitrm(data_table, 'placebo_pre-alprazolam_post ~ 1', 'WithinDesign', wd);

% check for approx. normal distribution
figure(1)
histogram(reshape(data_table{:,:},[],1), 10)

% plot as boxplot
figure(2)
boxplot(data_table{:,:})
set(gca, 'FontSize', 14)
set(gca, 'xtick', 1:4, 'xticklabel', {'pre' 'post' 'pre' 'post'})
yl = get(gca, 'ylim')
text(1.15, yl(1) - (yl(2) - yl(1))*0.1, 'placebo', 'FontSize', 14)
text(3.15, yl(1) - (yl(2) - yl(1))*0.1, 'alprazolam', 'FontSize', 14)

% check for sphericity
rm.mauchly

% results table
ranova_data = ranova(rm, 'WithinModel', 'medication*time');

% save ANOVA output 
ranova_data = ranova_data(:, 1:6);
writetable(ranova_data, [folder_output '\GABA_YC_ranova.xlsx'],  'Sheet', 'rsEEG - sigma')

% add mean values
writematrix(data_avg(:, :, 1), [folder_output '\GABA_YC_ranova.xlsx'],...
    'Sheet', 'rsEEG - sigma', 'WriteMode', 'append')
writematrix(data_avg(:, :, 3), [folder_output '\GABA_YC_ranova.xlsx'],...
    'Sheet', 'rsEEG - sigma', 'WriteMode', 'append')
clear data data_avg data_table rm wd ranova_data

%% 11) RANOVA: rsEEG - delta
% ----- extract data -----
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            statement = ['data(p, m, t) = GABA_YC_results.rsEEG(m).delta.' time{t} '(p);'];
            eval(statement)
        end
    end
end
clear m p t statement

% extract average values
for m = 1:length(medication)
    for t = 1:length(time)
        data_avg(m, t, 1) = mean(data(:, m, t), 'omitnan');
        data_avg(m, t, 2) = std(data(:, m, t), 'omitnan');      % STD
        n = sum(~isnan(data(:, m, t)));
        data_avg(m, t, 3) = data_avg(m, t, 2)/sqrt(n);              % SEM
        data_avg(m, t, 4) = data_avg(m, t, 3) * z;                  % 95% confidence interval      
    end
end
clear t m n
% pd = fitdist(data_med(1, :, 1)','Normal')
% CI = paramci(pd);
    
% ----- test for differences - rm-ANOVA -----
% prepare data table
data_table = table;
for m = 1:length(medication)
    for t = 1:length(time)
        statement = ['data_table.' medication{m} '_' time{t} ' = data(:, m, t);'];
        eval(statement)
    end
end
clear m t statement 

% within-subjects design
wd = table([repmat(medication(1), length(time), 1); repmat(medication(2), length(time), 1)],...
    [1:length(time) 1:length(time)]','VariableNames',{'medication','time'});
wd.medication = categorical(wd.medication); wd.time = categorical(wd.time);

% rm-ANOVA
rm = fitrm(data_table, 'placebo_pre-alprazolam_post ~ 1', 'WithinDesign', wd);

% check for approx. normal distribution
figure(1)
histogram(reshape(data_table{:,:},[],1), 10)

% plot as boxplot
figure(2)
boxplot(data_table{:,:})
set(gca, 'FontSize', 14)
set(gca, 'xtick', 1:4, 'xticklabel', {'pre' 'post' 'pre' 'post'})
yl = get(gca, 'ylim')
text(1.15, yl(1) - (yl(2) - yl(1))*0.1, 'placebo', 'FontSize', 14)
text(3.15, yl(1) - (yl(2) - yl(1))*0.1, 'alprazolam', 'FontSize', 14)

% check for sphericity
rm.mauchly

% results table
ranova_data = ranova(rm, 'WithinModel', 'medication*time');

% save ANOVA output 
ranova_data = ranova_data(:, 1:6);
writetable(ranova_data, [folder_output '\GABA_YC_ranova.xlsx'],  'Sheet', 'rsEEG - delta')

% add mean values
writematrix(data_avg(:, :, 1), [folder_output '\GABA_YC_ranova.xlsx'],...
    'Sheet', 'rsEEG - delta', 'WriteMode', 'append')
writematrix(data_avg(:, :, 3), [folder_output '\GABA_YC_ranova.xlsx'],...
    'Sheet', 'rsEEG - delta', 'WriteMode', 'append')
clear data data_avg data_table rm wd ranova_data

%% 12) RANOVA: rsEEG - AAC
% ----- extract data -----
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            statement = ['data(p, m, t) = GABA_YC_results.rsEEG(m).AAC.' time{t} '(p);'];
            eval(statement)
        end
    end
end
clear m p t statement

% extract average values
for m = 1:length(medication)
    for t = 1:length(time)
        data_avg(m, t, 1) = mean(data(:, m, t), 'omitnan');
        data_avg(m, t, 2) = std(data(:, m, t), 'omitnan');      % STD
        n = sum(~isnan(data(:, m, t)));
        data_avg(m, t, 3) = data_avg(m, t, 2)/sqrt(n);              % SEM
        data_avg(m, t, 4) = data_avg(m, t, 3) * z;                  % 95% confidence interval      
    end
end
clear t m n
% pd = fitdist(data_med(1, :, 1)','Normal')
% CI = paramci(pd);
    
% ----- test for differences - rm-ANOVA -----
% prepare data table
data_table = table;
for m = 1:length(medication)
    for t = 1:length(time)
        statement = ['data_table.' medication{m} '_' time{t} ' = data(:, m, t);'];
        eval(statement)
    end
end
clear m t statement 

% within-subjects design
wd = table([repmat(medication(1), length(time), 1); repmat(medication(2), length(time), 1)],...
    [1:length(time) 1:length(time)]','VariableNames',{'medication','time'});
wd.medication = categorical(wd.medication); wd.time = categorical(wd.time);

% rm-ANOVA
rm = fitrm(data_table, 'placebo_pre-alprazolam_post ~ 1', 'WithinDesign', wd);

% check for approx. normal distribution
figure(1)
histogram(reshape(data_table{:,:},[],1), 10)

% plot as boxplot
figure(2)
boxplot(data_table{:,:})
set(gca, 'FontSize', 14)
set(gca, 'xtick', 1:4, 'xticklabel', {'pre' 'post' 'pre' 'post'})
yl = get(gca, 'ylim')
text(1.15, yl(1) - (yl(2) - yl(1))*0.1, 'placebo', 'FontSize', 14)
text(3.15, yl(1) - (yl(2) - yl(1))*0.1, 'alprazolam', 'FontSize', 14)

% check for sphericity
rm.mauchly

% results table
ranova_data = ranova(rm, 'WithinModel', 'medication*time');

% save ANOVA output 
ranova_data = ranova_data(:, 1:6);
writetable(ranova_data, [folder_output '\GABA_YC_ranova.xlsx'],  'Sheet', 'rsEEG - AAC')

% add mean values
writematrix(data_avg(:, :, 1), [folder_output '\GABA_YC_ranova.xlsx'],...
    'Sheet', 'rsEEG - AAC', 'WriteMode', 'append')
writematrix(data_avg(:, :, 3), [folder_output '\GABA_YC_ranova.xlsx'],...
    'Sheet', 'rsEEG - AAC', 'WriteMode', 'append')
clear data data_avg data_table rm wd ranova_data yl

%% 13) RANOVA: rsEEG - SE
% ----- extract data -----
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            statement = ['data(p, m, t) = GABA_YC_results.rsEEG(m).SE.' time{t} '.closed(p);'];
            eval(statement)
        end
    end
end
clear m p t statement

% extract average values
for m = 1:length(medication)
    for t = 1:length(time)
        data_avg(m, t, 1) = mean(data(:, m, t), 'omitnan');
        data_avg(m, t, 2) = std(data(:, m, t), 'omitnan');      % STD
        n = sum(~isnan(data(:, m, t)));
        data_avg(m, t, 3) = data_avg(m, t, 2)/sqrt(n);              % SEM
        data_avg(m, t, 4) = data_avg(m, t, 3) * z;                  % 95% confidence interval      
    end
end
clear t m n
% pd = fitdist(data_med(1, :, 1)','Normal')
% CI = paramci(pd);
    
% ----- test for differences - rm-ANOVA -----
% prepare data table
data_table = table;
for m = 1:length(medication)
    for t = 1:length(time)
        statement = ['data_table.' medication{m} '_' time{t} ' = data(:, m, t);'];
        eval(statement)
    end
end
clear m t statement 

% within-subjects design
wd = table([repmat(medication(1), length(time), 1); repmat(medication(2), length(time), 1)],...
    [1:length(time) 1:length(time)]','VariableNames',{'medication','time'});
wd.medication = categorical(wd.medication); wd.time = categorical(wd.time);

% rm-ANOVA
rm = fitrm(data_table, 'placebo_pre-alprazolam_post ~ 1', 'WithinDesign', wd);

% check for approx. normal distribution
figure(1)
histogram(reshape(data_table{:,:},[],1), 10)

% plot as boxplot
figure(2)
boxplot(data_table{:,:})
set(gca, 'FontSize', 14)
set(gca, 'xtick', 1:4, 'xticklabel', {'pre' 'post' 'pre' 'post'})
yl = get(gca, 'ylim')
text(1.15, yl(1) - (yl(2) - yl(1))*0.1, 'placebo', 'FontSize', 14)
text(3.15, yl(1) - (yl(2) - yl(1))*0.1, 'alprazolam', 'FontSize', 14)

% check for sphericity
rm.mauchly

% results table
ranova_data = ranova(rm, 'WithinModel', 'medication*time');

% save ANOVA output 
ranova_data = ranova_data(:, 1:6);
writetable(ranova_data, [folder_output '\GABA_YC_ranova.xlsx'],  'Sheet', 'rsEEG - SE')

% add mean values
writematrix(data_avg(:, :, 1), [folder_output '\GABA_YC_ranova.xlsx'],...
    'Sheet', 'rsEEG - SE', 'WriteMode', 'append')
writematrix(data_avg(:, :, 3), [folder_output '\GABA_YC_ranova.xlsx'],...
    'Sheet', 'rsEEG - SE', 'WriteMode', 'append')
clear data data_avg data_table rm wd ranova_data yl

%% 14a) CORRELATION: TEP x MEP
% variable names
varnames = [peaks_M1, {'MEP'}];

% ----- extract data: all TS and ppTMS together -----
data_cor = [];
for m = 1:length(medication)
    for s = [2 3]
        % TEPs 
        for k = 1:length(peaks_M1)
            row_counter = 1;
            for p = 1:length(participant)               
                data_i(row_counter:row_counter+1, k) = [GABA_YC_results.TEP_M1(m).amplitude.pre(p, s, k); GABA_YC_results.TEP_M1(m).amplitude.post(p, s, k)]; 
                row_counter = row_counter + 2;
            end
        end

        % MEPs
        row_counter = 1;
        for p = 1:length(participant)              
            data_i(row_counter:row_counter+1, k+1) = [GABA_YC_results.MEP(m).amplitude.pre(p, s-1); GABA_YC_results.MEP(m).amplitude.post(p, s-1)];
            row_counter = row_counter + 2;
        end
        
        % append to global data matrix
        data_cor = [data_cor; data_i];
    end
end
clear m s p k data_i

% ----- overview of correlation between all variables - ranked data -----             
% create correlation matrix
[fig, coefficients, p_values] = correlation_preview(data_cor, varnames, 'method', 'Spearman'); 

% save figure
fig_name = 'corr_TEPxMEP_ranked_overview';
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name  '.png'])
figure_counter = figure_counter + 1;
clear fig fig_name

% ----- significant ranked correlation -----  
% Bonferroni correction of alpha
alpha_cor = alpha/length(peaks);

% identify significant cases        
[cor_coef, cor_p] = corr(data_cor, 'Type', 'Spearman');
[row, col] = find(cor_p < alpha_cor);

% rank the data
for a = 1:size(data_cor, 2)
    [temp, data_cor(:, a)]  = ismember(data_cor(:, a), unique(data_cor(:, a)));
end

% plot significant cases
for a = 1:length(row)    
    % prepare linear model: y ~ 1 + x
    data_model = fitlm(data_cor(:, row(a)), data_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

    % choose only correlations that show TEP-MEP interactions
    if ismember(row(a), 1:6) && col(a) == 7     
        % plot data + regression line
        fig = figure(figure_counter);
        hold on
        plot_cor = plotAdded(data_model);

        % adjust parameters    
        title(['Spearman rank correlation: ' varnames{col(a)} ' ~ ' varnames{row(a)}])
        xlabel(['change in ' varnames{row(a)}]); ylabel(['change in ' varnames{col(a)}]);
        set(gca, 'FontSize', 14)
        plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
        plot_cor(1).MarkerEdgeColor = colours(2, :); plot_cor(1).MarkerFaceColor = colours(2, :);
        plot_cor(2).Color = colours(3, :); plot_cor(2).LineWidth = 2; 
        plot_cor(3).Color = colours(3, :); plot_cor(3).LineWidth = 2;
        legend off
        if data_model.Coefficients.Estimate(2) > 0
            text_pos = [0.95 0.85 0.75];
        else
            text_pos = [0.25 0.15 0.05];
        end
        T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
        T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
        T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.5f', cor_coef(row(a), col(a)), cor_p(row(a), col(a))), 'Units', 'Normalized');
        set(T(1), 'fontsize', 14, 'fontweight', 'bold', 'fontangle', 'italic'); 
        set(T(2), 'fontsize', 14); 
        set(T(3), 'fontsize', 14, 'color', colours(3, :)); 

        % save figure and continue
%         fig_name = ['corr_TEPxMEP_ranked_' varnames{row(a)} '_' varnames{col(a)}];
%         savefig([folder_figures '\' fig_name '.fig'])
%         saveas(fig, [folder_figures '\' fig_name  '.svg'])  
%         figure_counter = figure_counter + 1;
    end
end
clear data_cor cor_coef cor_p data_model plot_cor row col fig figure_name T text_pos m s a temp row_counter 

% ----- extract data: only TS baseline -----
for m = 1:length(medication)
    % TEPs 
    for k = 1:length(peaks)
        for p = 1:length(participant)               
            data_cor((m-1)*length(participant) + p, k) = GABA_YC_results.TEP_M1(m).amplitude.pre(p, 2, k);
        end
    end

    % MEPs
    for p = 1:length(participant)               
        data_cor((m-1)*length(participant) + p, k + 1) = GABA_YC_results.MEP(m).amplitude.pre(p, 1);
    end
end
clear m p k 

% ----- overview of correlation between all variables - ranked data -----             
% create correlation matrix
fig = correlation_preview(data_cor, varnames, 'method', 'Pearson'); 
fig = correlation_preview(data_cor, varnames, 'method', 'Spearman'); 

% save figure
fig_name = 'corr_TEPxMEP_bl_ranked_overview';
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name  '.png'])   
figure_counter = figure_counter + 1;
clear varnames fig fig_name

%% 14b) CORRELATION: TEP x MEP - N17
varnames = {'N17' 'MEP'};

% get the data
data_cor = [];
for m = 1:length(medication)
    for s = [2 3]
        % TEPs 
        row_counter = 1;
        for p = 1:length(participant)               
            data_i(row_counter:row_counter+1, 1) = [GABA_YC_results.TEP_M1(m).amplitude.pre(p, s, 1); GABA_YC_results.TEP_M1(m).amplitude.post(p, s, 1)]; 
            row_counter = row_counter + 2;
        end

        % MEPs
        row_counter = 1;
        for p = 1:length(participant)              
            data_i(row_counter:row_counter+1, 2) = [GABA_YC_results.MEP(m).amplitude.pre(p, s-1); GABA_YC_results.MEP(m).amplitude.post(p, s-1)];
            row_counter = row_counter + 2;
        end
        
        % append to global data matrix
        data_cor = [data_cor; data_i];
    end
end
clear m s p k data_i

% % rank the data - reverse order of N17 --> more negative = bigger amplitude
% [sd, data_cor_ranked(:, 1)] = sort(data_cor(:, 1), 'descend');
% [sd, data_cor_ranked(:, 2)] = sort(data_cor(:, 2), 'ascend');
% clear sd

% rank the data
for a = 1:size(data_cor, 2)
    [temp, data_cor(:, a)]  = ismember(data_cor(:, a), unique(data_cor(:, a)));
end

% identify p    
[cor_coef, cor_p] = corr(data_cor, 'Type', 'Spearman');

% prepare linear model: y ~ 1 + x
data_model = fitlm(data_cor(:, 1), data_cor(:, 2), 'VarNames', varnames);

% plot data + regression line
fig = figure(figure_counter);
hold on
plot_cor = plotAdded(data_model);

% adjust parameters    
% title(['Spearman rank correlation: ' varnames{2} ' ~ ' varnames{1}])
xlabel(varnames{1}); ylabel(varnames{2});
set(gca, 'FontSize', 20)
plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
plot_cor(1).MarkerEdgeColor = colours(1, :); plot_cor(1).MarkerFaceColor = colours(1, :);
plot_cor(2).Color = [0 0 0]; plot_cor(2).LineWidth = 4; 
plot_cor(3).Color = [0 0 0]; plot_cor(3).LineWidth = 2;
legend off
text_pos = [0.90 0.78 0.66];
T(1) = text(0.1, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
T(2) = text(0.1, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
T(3) = text(0.1, text_pos(3), sprintf('p = %1.5f', cor_p(1, 2)), 'Units', 'Normalized');
set(T(1), 'fontsize', 20, 'fontangle', 'italic'); 
set(T(2), 'fontsize', 20, 'fontweight', 'bold'); 
set(T(3), 'fontsize', 20, 'fontweight', 'bold'); 
set(gca, 'XDir','reverse')
title ''

% save figure and continue
fig_name = 'corr_N17xMEP_all_ranked';
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name  '.svg'])  
figure_counter = figure_counter + 1;

%% 15a) CORRELATION: TEP change x RS-EEG change
% variable names
varnames = [peaks, {'beta' 'delta' 'AAC' 'SE'}];

% ----- extract data -----
for m = 1:length(medication)
    for s = [1 2]
        % change in TEPs 
        for k = 1:length(peaks)
            for p = 1:length(participant)               
                data_cor(m, s, p, k) = GABA_YC_results.TEP_M1(m).amplitude.change(p, s, k);           
            end
        end


        % change in rsEEG measures
        for p = 1:length(participant)              
            data_cor(m, s, p, k+1) = GABA_YC_results.rsEEG(m).sigma.change(p);   
            data_cor(m, s, p, k+2) = GABA_YC_results.rsEEG(m).delta.change(p);
            data_cor(m, s, p, k+3) = GABA_YC_results.rsEEG(m).AAC.change(p);
            data_cor(m, s, p, k+4) = GABA_YC_results.rsEEG(m).SE.change.closed(p);
        end
    end
end
clear m s p k

% ----- overview of linear correlation between all variables -----
for m = 1:length(medication)
    for s = 2%[1 2]
        % select the data                
        mat_cor = squeeze(data_cor(m, s, :, :));
        
        % create correlation matrix
        fig = correlation_preview(mat_cor, varnames, 'method', 'Pearson'); 
        
        % save figure
%         fig_name = ['corr_TEPxrsEEG_linear_' medication{m} '_' stimulus{s} '_overview'];
%         savefig([folder_figures '\' fig_name '.fig'])
%         saveas(fig, [folder_figures '\' fig_name  '.png'])   
%         
        % update figure counter
        figure_counter = figure_counter + 1;
    end
end
clear m s fig fig_name

% ----- overview of correlation between all variables - ranked data -----
for m = 1:length(medication)
    for s = 2%[1 2]
        % select the data                
        mat_cor = squeeze(data_cor(m, s, :, :));
        
        % create correlation matrix
        fig = correlation_preview(mat_cor, varnames, 'method', 'Spearman'); 
        
        % save figure
%         fig_name = ['corr_TEPxrsEEG_ranked_' medication{m} '_' stimulus{s} '_overview'];
%         savefig([folder_figures '\' fig_name '.fig'])
%         saveas(fig, [folder_figures '\' fig_name  '.png']) 
        
        % update figure counter
        figure_counter = figure_counter + 1;
    end
end
clear m s fig fig_name

% ----- significant linear correlation -----
for m = 1:length(medication)
    for s = [1 2]
        % Bonferroni correction of alpha
        if s == 1
            peaks_n = 5;
        else
            peaks_n = 6;
        end
        alpha_cor = alpha/(peaks_n * 4);
        
        % select the data                
        mat_cor = squeeze(data_cor(m, s, :, :));
        
        % identify significant cases        
        [cor_coef, cor_p] = corrcoef(mat_cor);
        [row, col] = find(cor_p < alpha_cor);
        
        % plot significant cases
        for a = 1:length(row)    
            % prepare linear model: y ~ 1 + x
            data_model = fitlm(mat_cor(:, row(a)), mat_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);
            
            % choose only correlations that show TEP-rsEEG interactions
            if ismember(row(a), 1:6) & ismember(col(a), 7:10)            
                % plot data + regression line
                fig = figure(figure_counter);
                hold on
                plot_cor = plotAdded(data_model);

                % adjust parameters    
                title(sprintf('Linear correlation - %s, %s\n%s ~ %s', medication{m}, stimulus{s}, varnames{col(a)}, varnames{row(a)}))
                xlabel(['change in ' varnames{row(a)}]); ylabel(['change in ' varnames{col(a)}]);
                set(gca, 'FontSize', 14)
                plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
                plot_cor(1).MarkerEdgeColor = colours(2, :); plot_cor(1).MarkerFaceColor = colours(2, :);
                plot_cor(2).Color = colours(3, :); plot_cor(2).LineWidth = 2; 
                plot_cor(3).Color = colours(3, :); plot_cor(3).LineWidth = 2;
                legend off
                if data_model.Coefficients.Estimate(2) > 0
                    text_pos = [0.95 0.85 0.75];
                else
                    text_pos = [0.25 0.15 0.05];
                end
                T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
                T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
                T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.5f', cor_coef(row(a), col(a)), cor_p(row(a), col(a))), 'Units', 'Normalized');
                set(T(1), 'fontsize', 14, 'fontweight', 'bold', 'fontangle', 'italic'); 
                set(T(2), 'fontsize', 14); 
                set(T(3), 'fontsize', 14, 'color', colours(3, :)); 

                % save figure and continue
                fig_name = ['corr_TEPxrsEEG_linear_' medication{m} '_' stimulus{s} '_' varnames{row(a)} '_' varnames{col(a)}];
                savefig([folder_figures '\' fig_name '.fig'])
                saveas(fig, [folder_figures '\' fig_name  '.png'])  
                figure_counter = figure_counter + 1;
            end
        end
    end
end
clear mat_cor cor_coef cor_p data_model plot_cor row col fig fig_name T text_pos m s a temp peaks_n alpha_cor

% ----- significant ranked correlation -----
for m = 1:length(medication)
    for s = 2     
        % Bonferroni correction of alpha
        if s == 1
            peaks_n = 5;
        else
            peaks_n = 6;
        end
        alpha_cor = alpha/(peaks_n * 4);
        
        % select the data                
        mat_cor = squeeze(data_cor(m, s, :, :));
        
        % identify significant cases        
        [cor_coef, cor_p] = corr(mat_cor, 'Type', 'Spearman');
        [row, col] = find(cor_p < alpha_cor);
        
        % rank the data
        for a = 1:size(mat_cor, 2)
            [temp, mat_cor(:, a)]  = ismember(mat_cor(:, a), unique(mat_cor(:, a)));
        end
        
        % plot significant cases
        for a = 1:length(row)    
            % prepare linear model: y ~ 1 + x
            data_model = fitlm(mat_cor(:, row(a)), mat_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);
            
            % choose only correlations that show TEP-MEP interactions
            if ismember(row(a), 1:6) & ismember(col(a), 7:10)           
                % plot data + regression line
                fig = figure(figure_counter);
                hold on
                plot_cor = plotAdded(data_model);

                % adjust parameters    
                title(sprintf('Spearman rank correlation - %s, %s\n%s ~ %s', medication{m}, stimulus{s}, varnames{col(a)}, varnames{row(a)}))
                xlabel(['change in ' varnames{row(a)}]); ylabel(['change in ' varnames{col(a)}]);
                set(gca, 'FontSize', 14)
                plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
                plot_cor(1).MarkerEdgeColor = colours(2, :); plot_cor(1).MarkerFaceColor = colours(2, :);
                plot_cor(2).Color = colours(3, :); plot_cor(2).LineWidth = 2; 
                plot_cor(3).Color = colours(3, :); plot_cor(3).LineWidth = 2;
                legend off
                if data_model.Coefficients.Estimate(2) > 0
                    text_pos = [0.95 0.85 0.75];
                else
                    text_pos = [0.25 0.15 0.05];
                end
                T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
                T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
                T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.5f', cor_coef(row(a), col(a)), cor_p(row(a), col(a))), 'Units', 'Normalized');
                set(T(1), 'fontsize', 14, 'fontweight', 'bold', 'fontangle', 'italic'); 
                set(T(2), 'fontsize', 14); 
                set(T(3), 'fontsize', 14, 'color', colours(3, :)); 

                % save figure and continue
                fig_name = ['corr_TEPxrsEEG_ranked_' medication{m} '_' stimulus{s} '_' varnames{row(a)} '_' varnames{col(a)}];
                savefig([folder_figures '\' fig_name '.fig'])
                saveas(fig, [folder_figures '\' fig_name  '.png'])  
                figure_counter = figure_counter + 1;
            end
        end
    end
end
clear varnames data_cor mat_cor cor_coef cor_p data_model plot_cor row col fig fig_name T text_pos m s a temp peaks_n alpha_cor

%% 16) CORRELATION: TEP SICI x MEP SICI 
% variable names
varnames = {'TEP-SICI' 'MEP-SICI'};

% ----- extract data: baseline -----
for m = 1:length(medication)
    for p = 1:length(participant) 
        % TEP-SICI 
        data_cor((m-1)*length(participant) + p, 1) = GABA_YC_results.SICI(m).GFP_AUC.pre(p);

        % MEP-SICI
        data_cor((m-1)*length(participant) + p, 2) = GABA_YC_results.SICI(m).MEP.pre(p);
    end
end
clear m p 

% ----- significant linear correlation -----
% identify significant cases        
[cor_coef, cor_p] = corrcoef(data_cor);
[row, col] = find(cor_p < alpha);

% plot significant cases
for a = 1:length(row)/2    
    % prepare linear model: y ~ 1 + x
    data_model = fitlm(data_cor(:, row(a)), data_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);
          
    % plot data + regression line
    fig = figure(figure_counter);
    hold on
    plot_cor = plotAdded(data_model);

    % adjust parameters    
    title(sprintf('Linear correlation - baseline:\n%s ~ %s', varnames{col(a)}, varnames{row(a)}))
    xlabel(['change in ' varnames{row(a)}]); ylabel(['change in ' varnames{col(a)}]);
    set(gca, 'FontSize', 14)
    plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
    plot_cor(1).MarkerEdgeColor = colours(2, :); plot_cor(1).MarkerFaceColor = colours(2, :);
    plot_cor(2).Color = colours(3, :); plot_cor(2).LineWidth = 2; 
    plot_cor(3).Color = colours(3, :); plot_cor(3).LineWidth = 2;
    legend off
    if data_model.Coefficients.Estimate(2) > 0
        text_pos = [0.95 0.85 0.75];
    else
        text_pos = [0.25 0.15 0.05];
    end
    T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
    T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
    T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.5f', cor_coef(row(a), col(a)), cor_p(row(a), col(a))), 'Units', 'Normalized');
    set(T(1), 'fontsize', 14, 'fontweight', 'bold', 'fontangle', 'italic'); 
    set(T(2), 'fontsize', 14); 
    set(T(3), 'fontsize', 14, 'color', colours(3, :)); 

    % save figure and continue
    fig_name = ['corr_TEP-SICI_GFP_AUCxMEP-SICI_linear_' varnames{row(a)} '_' varnames{col(a)}];
    savefig([folder_figures '\' fig_name '.fig'])
    saveas(fig, [folder_figures '\' fig_name  '.png'])  
    figure_counter = figure_counter + 1;
end
clear cor_coef alpha_cor cor_p data_model plot_cor row col fig fig_name T text_pos a temp alpha_cor

% ----- significant ranked correlation -----   
% identify significant cases        
[cor_coef, cor_p] = corr(data_cor, 'Type', 'Spearman');
[row, col] = find(cor_p < alpha);

% rank the data
for a = 1:size(data_cor, 2)
    [temp, data_cor(:, a)]  = ismember(data_cor(:, a), unique(data_cor(:, a)));
end

% plot significant cases
for a = 1:length(row)/2    
    % prepare linear model: y ~ 1 + x
    data_model = fitlm(data_cor(:, row(a)), data_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

    % plot data + regression line
    fig = figure(figure_counter);
    hold on
    plot_cor = plotAdded(data_model);

    % adjust parameters    
    title(sprintf('Spearman rank correlation - baseline:\n%s ~ %s', varnames{col(a)}, varnames{row(a)}))
    xlabel(['change in ' varnames{row(a)}]); ylabel(['change in ' varnames{col(a)}]);
    set(gca, 'FontSize', 14)
    plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
    plot_cor(1).MarkerEdgeColor = colours(2, :); plot_cor(1).MarkerFaceColor = colours(2, :);
    plot_cor(2).Color = colours(3, :); plot_cor(2).LineWidth = 2; 
    plot_cor(3).Color = colours(3, :); plot_cor(3).LineWidth = 2;
    legend off
    if data_model.Coefficients.Estimate(2) > 0
        text_pos = [0.95 0.85 0.75];
    else
        text_pos = [0.25 0.15 0.05];
    end
    T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
    T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
    T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.5f', cor_coef(row(a), col(a)), cor_p(row(a), col(a))), 'Units', 'Normalized');
    set(T(1), 'fontsize', 14, 'fontweight', 'bold', 'fontangle', 'italic'); 
    set(T(2), 'fontsize', 14); 
    set(T(3), 'fontsize', 14, 'color', colours(3, :)); 

    % save figure and continue
    fig_name = ['corr_TEP-SICI_GFP_AUCxMEP-SICI_ranked_' varnames{row(a)} '_' varnames{col(a)}];
    savefig([folder_figures '\' fig_name '.fig'])
    saveas(fig, [folder_figures '\' fig_name  '.png'])  
    figure_counter = figure_counter + 1;
end
clear varnames data_cor cor_coef cor_p data_model plot_cor row col fig fig_name T text_pos m s a temp peaks_n alpha_cor

%% 17) CORRELATION: TEP SICI x MEP SICI 
% variable names
varnames = {'TEP-SICI P1' 'TEP-SICI N2' 'MEP-SICI'};

% ----- extract data: baseline -----
for m = 1:length(medication)
    % TEP-SICI 
    for k = [1 2]
        for p = 1:length(participant)               
            data_cor((m-1)*length(participant) + p, k) = GABA_YC_results.SICI(m).TEP_SICI.pre(p, k);
        end
    end

    % MEP-SICI
    for p = 1:length(participant)               
        data_cor((m-1)*length(participant) + p, k + 1) = GABA_YC_results.SICI(m).MEP.pre(p);
    end
end
clear m p k 

% ----- significant linear correlation -----
% Bonferroni correction of alpha
alpha_cor = alpha/2;

% identify significant cases        
[cor_coef, cor_p] = corrcoef(data_cor);
[row, col] = find(cor_p < alpha_cor);

% plot significant cases
for a = 1:length(row)    
    % prepare linear model: y ~ 1 + x
    data_model = fitlm(data_cor(:, row(a)), data_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

    % choose only correlations that show TEP-rsEEG interactions
    if ismember(row(a), [1 2]) && col(a) == 3            
        % plot data + regression line
        fig = figure(figure_counter);
        hold on
        plot_cor = plotAdded(data_model);

        % adjust parameters    
        title(sprintf('Linear correlation - baseline:\n%s ~ %s', varnames{col(a)}, varnames{row(a)}))
        xlabel(['change in ' varnames{row(a)}]); ylabel(['change in ' varnames{col(a)}]);
        set(gca, 'FontSize', 14)
        plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
        plot_cor(1).MarkerEdgeColor = colours(2, :); plot_cor(1).MarkerFaceColor = colours(2, :);
        plot_cor(2).Color = colours(3, :); plot_cor(2).LineWidth = 2; 
        plot_cor(3).Color = colours(3, :); plot_cor(3).LineWidth = 2;
        legend off
        if data_model.Coefficients.Estimate(2) > 0
            text_pos = [0.95 0.85 0.75];
        else
            text_pos = [0.25 0.15 0.05];
        end
        T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
        T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
        T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.5f', cor_coef(row(a), col(a)), cor_p(row(a), col(a))), 'Units', 'Normalized');
        set(T(1), 'fontsize', 14, 'fontweight', 'bold', 'fontangle', 'italic'); 
        set(T(2), 'fontsize', 14); 
        set(T(3), 'fontsize', 14, 'color', colours(3, :)); 

        % save figure and continue
        fig_name = ['corr_TEP-SICIxMEP-SICI_linear_' varnames{row(a)} '_' varnames{col(a)}];
        savefig([folder_figures '\' fig_name '.fig'])
        saveas(fig, [folder_figures '\' fig_name  '.png'])  
        figure_counter = figure_counter + 1;
    end
end
clear cor_coef alpha_cor cor_p data_model plot_cor row col fig fig_name T text_pos a temp alpha_cor

% ----- significant ranked correlation -----   
% Bonferroni correction of alpha
alpha_cor = alpha/2;

% identify significant cases        
[cor_coef, cor_p] = corr(data_cor, 'Type', 'Spearman');
[row, col] = find(cor_p < alpha_cor);

% rank the data
for a = 1:size(data_cor, 2)
    [temp, data_cor(:, a)]  = ismember(data_cor(:, a), unique(data_cor(:, a)));
end

% plot significant cases
for a = 1:length(row)    
    % prepare linear model: y ~ 1 + x
    data_model = fitlm(data_cor(:, row(a)), data_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

    % choose only correlations that show TEP-MEP interactions
    if ismember(row(a), [1 2]) && col(a) == 3           
        % plot data + regression line
        fig = figure(figure_counter);
        hold on
        plot_cor = plotAdded(data_model);

        % adjust parameters    
        title(sprintf('Spearman rank correlation - baseline:\n%s ~ %s', varnames{col(a)}, varnames{row(a)}))
        xlabel(['change in ' varnames{row(a)}]); ylabel(['change in ' varnames{col(a)}]);
        set(gca, 'FontSize', 14)
        plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
        plot_cor(1).MarkerEdgeColor = colours(2, :); plot_cor(1).MarkerFaceColor = colours(2, :);
        plot_cor(2).Color = colours(3, :); plot_cor(2).LineWidth = 2; 
        plot_cor(3).Color = colours(3, :); plot_cor(3).LineWidth = 2;
        legend off
        if data_model.Coefficients.Estimate(2) > 0
            text_pos = [0.95 0.85 0.75];
        else
            text_pos = [0.25 0.15 0.05];
        end
        T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
        T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
        T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.5f', cor_coef(row(a), col(a)), cor_p(row(a), col(a))), 'Units', 'Normalized');
        set(T(1), 'fontsize', 14, 'fontweight', 'bold', 'fontangle', 'italic'); 
        set(T(2), 'fontsize', 14); 
        set(T(3), 'fontsize', 14, 'color', colours(3, :)); 

        % save figure and continue
        fig_name = ['corr_TEP-SICIxMEP-SICI_ranked_' varnames{row(a)} '_' varnames{col(a)}];
        savefig([folder_figures '\' fig_name '.fig'])
        saveas(fig, [folder_figures '\' fig_name  '.png'])  
        figure_counter = figure_counter + 1;
    end
end
clear varnames data_cor cor_coef cor_p data_model plot_cor row col fig fig_name T text_pos m s a temp peaks_n alpha_cor

%% 18) CORRELATION: TEP-peak SICI x MEP SICI 
% variable names
varnames = [peaks_M1 {'MEP-SICI'}];

% ----- extract data: baseline -----
for m = 1:length(medication)
    % TEP peak-SICI 
    for k = 1:length(peaks_M1)
        for p = 1:length(participant)               
            data_cor((m-1)*length(participant) + p, k) = GABA_YC_results.SICI(m).TEP.pre(p, k);
        end
    end

    % MEP-SICI
    for p = 1:length(participant)               
        data_cor((m-1)*length(participant) + p, k + 1) = GABA_YC_results.SICI(m).MEP.pre(p);
    end
end
clear m p k 

% ----- significant linear correlation -----
% Bonferroni correction of alpha
alpha_cor = alpha/length(peaks_M1);

% identify significant cases        
[cor_coef, cor_p] = corrcoef(data_cor);
[row, col] = find(cor_p < alpha_cor);

% plot significant cases
for a = 1:length(row)    
    % prepare linear model: y ~ 1 + x
    data_model = fitlm(data_cor(:, row(a)), data_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

    % choose only correlations that show TEP-rsEEG interactions
    if ismember(row(a), 1:length(peaks_M1)) && col(a) == length(peaks_M1) + 1            
        % plot data + regression line
        fig = figure(figure_counter);
        hold on
        plot_cor = plotAdded(data_model);

        % adjust parameters    
        title(sprintf('Linear correlation - baseline:\n%s ~ %s', varnames{col(a)}, varnames{row(a)}))
        xlabel(['change in ' varnames{row(a)}]); ylabel(varnames{col(a)});
        set(gca, 'FontSize', 14)
        plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
        plot_cor(1).MarkerEdgeColor = colours(2, :); plot_cor(1).MarkerFaceColor = colours(2, :);
        plot_cor(2).Color = colours(3, :); plot_cor(2).LineWidth = 2; 
        plot_cor(3).Color = colours(3, :); plot_cor(3).LineWidth = 2;
        legend off
        if data_model.Coefficients.Estimate(2) > 0
            text_pos = [0.95 0.85 0.75];
        else
            text_pos = [0.25 0.15 0.05];
        end
        T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
        T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
        T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.5f', cor_coef(row(a), col(a)), cor_p(row(a), col(a))), 'Units', 'Normalized');
        set(T(1), 'fontsize', 14, 'fontweight', 'bold', 'fontangle', 'italic'); 
        set(T(2), 'fontsize', 14); 
        set(T(3), 'fontsize', 14, 'color', colours(3, :)); 

        % save figure and continue
        fig_name = ['corr_TEP-SICIxMEP-SICI_linear_' varnames{row(a)} '_' varnames{col(a)}];
        savefig([folder_figures '\' fig_name '.fig'])
        saveas(fig, [folder_figures '\' fig_name  '.png'])  
        figure_counter = figure_counter + 1;
    end
end
clear cor_coef alpha_cor cor_p data_model plot_cor row col fig fig_name T text_pos a temp alpha_cor 

% ----- significant ranked correlation -----  
% Bonferroni correction of alpha
alpha_cor = alpha/length(peaks_M1);

% identify significant cases        
[cor_coef, cor_p] = corr(data_cor, 'Type', 'Spearman');
[row, col] = find(cor_p < alpha_cor);

% rank the data
for a = 1:size(data_cor, 2)
    [temp, data_cor(:, a)]  = ismember(data_cor(:, a), unique(data_cor(:, a)));
end

% plot significant cases
for a = 1:length(row)    
    % prepare linear model: y ~ 1 + x
    data_model = fitlm(data_cor(:, row(a)), data_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

    % choose only correlations that show TEP-rsEEG interactions
    if ismember(row(a), 1:length(peaks_M1)) && col(a) == length(peaks_M1) + 1            
        % plot data + regression line
        fig = figure(figure_counter);
        hold on
        plot_cor = plotAdded(data_model);

        % adjust parameters    
        title(sprintf('Spearman rank correlation - baseline:\n%s ~ %s', varnames{col(a)}, varnames{row(a)}))
        xlabel(['change in ' varnames{row(a)}]); ylabel(varnames{col(a)});
        set(gca, 'FontSize', 14)
        plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
        plot_cor(1).MarkerEdgeColor = colours(2, :); plot_cor(1).MarkerFaceColor = colours(2, :);
        plot_cor(2).Color = colours(3, :); plot_cor(2).LineWidth = 2; 
        plot_cor(3).Color = colours(3, :); plot_cor(3).LineWidth = 2;
        legend off
        if data_model.Coefficients.Estimate(2) > 0
            text_pos = [0.95 0.85 0.75];
        else
            text_pos = [0.25 0.15 0.05];
        end
        T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
        T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
        T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.5f', cor_coef(row(a), col(a)), cor_p(row(a), col(a))), 'Units', 'Normalized');
        set(T(1), 'fontsize', 14, 'fontweight', 'bold', 'fontangle', 'italic'); 
        set(T(2), 'fontsize', 14); 
        set(T(3), 'fontsize', 14, 'color', colours(3, :)); 

        % save figure and continue
        fig_name = ['corr_TEP-SICIxMEP-SICI_ranked_' varnames{row(a)} '_' varnames{col(a)}];
        savefig([folder_figures '\' fig_name '.fig'])
        saveas(fig, [folder_figures '\' fig_name  '.png'])  
        figure_counter = figure_counter + 1;
    end
end
clear varnames data_cor cor_coef alpha_cor cor_p data_model plot_cor row col fig fig_name T text_pos a temp alpha_cor 

