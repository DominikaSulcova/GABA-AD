%% GABA-AD: GROUP STATISTICS
% Written by Dominika for GABA-AD project (2021)
% 
% ----- DATA EXTRACTION -----
% Extracts subject info, parameters and outcome measures from all available 
% datasets and saves it in structure 'GABA_results'
% 1) subject info
%       - prepares outcome structure: one row per participant
%       - fills in fields 'subject' and 'session' --> rMT, vigilance etc.
% 2) TEP
%       - dataset info --> number of ICs removed during ICA
%       - mean peak amplitudes and latencies
%       - change in peak amplitude 
%       - SICI
% 3) RS-EEG
%       - dataset info --> ICs, epochs etc.
%       - data parameters --> IAF, TF, fbands 
%       - mean amplitude over fbands + amplitude change
%       - alpha attenuation coeficient
%       - spectral exponent 
%       - broad beta amplitude over central region + change    
% 4) MEP
%       - peak-to-peak amplitude 
%       - amplitude change
% 5) calculate mean values
% 6) export data for R
% 
% ----- DESCRIPTIVE STATISTICS -----
% 7) effect of medication 
% 
% ----- CORRELATIONS -----
% ) TEP x RS-EEG correlation
%       - calculates correlation coefficients and p-values for following variables:
%           - change in amplitude for all TEP peaks 
%           - change in AAC - only high alpha 3
%           - change in SE (eyes open and closed) - only 0.1 - 20Hz
%       - plots a correlation matrix for each medication and stimulus, 
%         visualize variable pairs that are significantly correlated 
%           --> uses fcn correlation_preview.m
%       - plots significantly correlated variables, shows regression line
%       - checks both linear and non-linear correlation (--> Spearman coefficient)
% 
% ) TEP x MEP correlation
%       - calculates correlation coefficients and p-values for following variables:
%           - change in amplitude for all TEP peaks 
%           - change in MEP peak-to-peak amplitude
%       - identifies significant correlations between MEP and any TEP peak change 
%       - plots significantly correlated variables, shows regression line
%       - checks both linear and non-linear correlation (--> Spearman coefficient)
% 
% ) TEP SICI x MEP SICI correlation
%       - only takes into an account data from the baseline
%       - calculates correlation coefficients and p-values for following variables:
%           - baseline SICI for all TEP peaks 
%           - baseline SICI in MEPs 
%       - identifies significant correlations between MEP and any TEP SICI
%       - plots significantly correlated variables, shows regression line
%       - checks both linear and non-linear correlation (--> Spearman coefficient)
% 
% ) TEP SICI x TEP change correlation
%       - calculates correlation coefficients and p-values for following variables:
%           - baseline SICI for all TEP peaks 
%           - change in TS TEP amplitude caused by the medication 
%       - identifies significant correlations between baseline TEP SICI and
%       alprazolam-induced change in the same TEP peak
%       - plots significantly correlated variables, shows regression line
%       - checks both linear and non-linear correlation (--> Spearman coefficient)


%% parameters
clear all, clc

% dataset
participant = 1:20;
stimulus = {'CS' 'TS' 'ppTMS'};
medication = {'placebo' 'alprazolam'};
time = {'pre' 'post'};

% stats, graphics
peaks = {'P30' 'N45' 'P60' 'N100' 'P180'};
alpha = 0.05;
z = 1.96;
load('colours2.mat')
shade = 0.2;
figure_counter = 1;

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
params = readtable('YC_parameters.csv');

% fill in basic info
if ~exist('GABA_results.mat')
    GABA_results = struct;
    for p = 1:length(participant)
        % subject info
        GABA_results(p).subject.age = params{p, 'Age'}; 
        GABA_results(p).subject.sex = params{p, 'Sex'};                         % 1 = female; 0 = male
        GABA_results(p).subject.handedness = params{p, 'Handedness'};           % 1 = right; 0 = left
        GABA_results(p).subject.APOE = params{p, 'APOE'};                       % 1 = homo/hetero eta4; 0 = no eta4

        % session info
        GABA_results(p).session(1).medication = 'placebo';
        GABA_results(p).session(1).rmt.pre = params{p, 'zyrtec_rMT_pre'}; 
        GABA_results(p).session(1).rmt.post = params{p, 'zyrtec_rMT_post'}; 
        GABA_results(p).session(1).rmt.change = params{p, 'zyrtec_rMT_post'} - params{p, 'xanax_rMT_pre'}; 
        GABA_results(p).session(1).vigilance(a) = params{p, 14 + a};

        GABA_results(p).session(2).medication = 'alprazolam';
        GABA_results(p).session(2).rmt.pre = params{p, 'xanax_rMT_pre'}; 
        GABA_results(p).session(2).rmt.post = params{p, 'xanax_rMT_post'}; 
        GABA_results(p).session(2).rmt.change = params{p, 'xanax_rMT_post'} - params{p, 'xanax_rMT_pre'}; 
        GABA_results(p).session(2).vigilance(a) = params{p, 7 + a};
    end
else
    load('GABA_results.mat')
end

%% 2) extract TEP measures
% fill in TEP related info
for p = 1:length(participant)
%     GABA_results(p).TEP(1).medication = 'placebo';
%     GABA_results(p).TEP(1).ICA.pre = params{p, 'zyrtec_ICA_pre'};
%     GABA_results(p).TEP(1).ICA.post = params{p, 'zyrtec_ICA_post'};
      GABA_results(p).TEP(1).epochs.pre(:) = squeeze(epochs(p, 1, 1, :))';
      GABA_results(p).TEP(1).epochs.post(:) = squeeze(epochs(p, 1, 2, :))';
      
%     GABA_results(p).TEP(2).medication = 'alprazolam';
%     GABA_results(p).TEP(2).ICA.pre = params{p, 'xanax_ICA_pre'};
%     GABA_results(p).TEP(2).ICA.post = params{p, 'xanax_ICA_post'};
      GABA_results(p).TEP(2).epochs.pre(:) = squeeze(epochs(p, 2, 1, :))';
      GABA_results(p).TEP(2).epochs.post(:) = squeeze(epochs(p, 2, 2, :))';
end
clear params epochs

% mean peak amplitudes and latencies
load('TEPs_final.mat')
for p = 1:length(participant)
    for m = 1:length(medication)
        for s = 1:length(stimulus)
            for t = 1:length(time)
                % choose data --> peak amplitudes and latencies
                rows = (TEPs_final.subject == p & ...
                    categorical(TEPs_final.medication) == medication{m} & ...
                    categorical(TEPs_final.stimulus) == stimulus{s} & ...
                    categorical(TEPs_final.time) == time{t});
                amps = TEPs_final.amplitude(rows)';
                lats = TEPs_final.latency(rows)';
                
                % append to the structure
                statement = ['GABA_results(p).TEP(m).' stimulus{s} '.' time{t} '(1, :) = amps;'];
                eval(statement)
                statement = ['GABA_results(p).TEP(m).' stimulus{s} '.' time{t} '(2, :) = lats;'];
                eval(statement)
            end
        end
    end
end
clear rows amps lats statement TEPs_final

% amplitude change
TEP_change = readtable('TEP_change_long.csv');
for p = 1:length(participant)
    for m = 1:length(medication)
        for s = 1:length(stimulus)
            for t = 1:length(time)
                % choose data 
                rows = (TEP_change.subject == p & ...
                    categorical(TEP_change.medication) == medication{m} & ...
                    categorical(TEP_change.stimulus) == stimulus{s});
                amp_change = TEP_change.amplitude(rows)';
                
                % append to the structure
                statement = ['GABA_results(p).TEP(m).' stimulus{s} '.change(1, :) = amp_change;'];
                eval(statement)
            end
        end
    end
end
clear rows amp_change statement TEP_change

% SICI
SICI = readtable('SICI_long.csv');
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            % choose data --> peak amplitudes and latencies
            rows = (SICI.subject == p & ...
                categorical(SICI.medication) == medication{m} & ...
                categorical(SICI.time) == time{t});
            sici_i = SICI.amplitude(rows)';

            % append to the structure
            statement = ['GABA_results(p).TEP(m).SICI.' time{t} '(1, :) = sici_i;'];
            eval(statement)
        end
    end
end
clear rows sici_i statement SICI

% save
save('GABA_results.mat', 'GABA_results')

%% 3) extract RS-EEG measures
% dataset info 
for p = 1:length(participant)
    for m = 1:length(medication)
        GABA_results(p).rsEEG(m).medication = medication{m};
        GABA_results(p).rsEEG(m).ROI = ROI;
    end
end

% IAF, TF, fband amplitude 
load('IAF.mat'); load('TF.mat'); load('fband.mat'); load('rsEEG_fband_amplitude.mat')
for p = 1:length(participant)
    for m = 1:length(medication)
        % fband names
        GABA_results(p).rsEEG(m).fband.bands = fband_names; 
        
        % choose data --> fband limits
        rows = (fband.subject == p & ...
            categorical(fband.medication) == medication{m});
        lim_i = fband{rows, 3:10};

        % append
        statement = ['GABA_results(p).rsEEG(m).fband.limits = lim_i;'];
        eval(statement)   
        
        % choose data --> TF
        rows = (TF.subject == p & ...
            categorical(TF.medication) == medication{m});
        tf_i = TF.occipital(rows);

        % append
        GABA_results(p).rsEEG(m).TF = tf_i;        
        
        for t = 1:length(time)
            % choose data --> IAF
            rows = (IAF.subject == p & ...
                categorical(IAF.medication) == medication{m} & ...
                categorical(IAF.time) == time{t});
            iaf_i = IAF{rows, 4:8};
            
            % append
            statement = ['GABA_results(p).rsEEG(m).IAF.' time{t} '(1, :) = iaf_i;'];
            eval(statement)      
            
            for c = 1:length(condition)
                for f = 1:length(fband_names)
                    % choose data --> fband amplitude
                    rows = (fband_amplitude.subject == p & ...
                        categorical(fband_amplitude.medication) == medication{m} & ...
                        categorical(fband_amplitude.time) == time{t} & ...
                        categorical(fband_amplitude.condition) == condition{c} & ...
                        categorical(fband_amplitude.fband) == fband_names{f});
                    amp_i = fband_amplitude{rows, 6:10}';

                    % append
                    statement = ['GABA_results(p).rsEEG(m).fband.amplitude.' time{t} '.' condition{c} '(:, f) = amp_i;'];
                    eval(statement)
                
                end
            end 
        end
        
        for c = 1:length(condition)
            for f = 1:length(fband_names)
                % choose data --> amplitude change
                for t = 1:length(time)
                    statement = ['data_' time{t} ' = GABA_results(p).rsEEG(m).fband.amplitude.' time{t} '.' condition{c} '(:, f);'];
                    eval(statement)
                end
                change_i = data_post - data_pre;

                % append
                statement = ['GABA_results(p).rsEEG(m).fband.amplitude.change.' condition{c} '(:, f) = change_i;'];
                eval(statement)

            end
        end         
    end
end
clear IAF iaf_i TF tf_i fband lim_i fband_amplitude amp_i change_i data_pre data_post rows statement

% alpha attenuation coeficient
load('AAC.mat'); load('AAC_change.mat');
for p = 1:length(participant)
    for m = 1:length(medication)
        GABA_results(p).rsEEG(m).AAC.bands = {'low alpha 1' 'low alpha 2'  'high alpha 3' 'broad alpha'};
        for t = 1:length(time)
            % AAC
            statement = ['GABA_results(p).rsEEG(m).AAC.' time{t} '(1, :) = squeeze(ACC(m, t, :, p));'];
            eval(statement)            
        end
        
        % AAC change
        GABA_results(p).rsEEG(m).AAC.change(1, :) = squeeze(ACC_change(m, :, p));     
    end
end
clear ACC ACC_change statement

% spectral exponent  
load('rsEEG_spect_exp.mat');
for p = 1:length(participant)
    for m = 1:length(medication)
        GABA_results(p).rsEEG(m).SE.bands = {'low' 'high' 'broad'};        
        for t = 1:length(time)
            for c = 1:length(condition)
                % choose data --> SE
                for a = 1:numel(GABA_results(p).rsEEG(m).SE.bands)
                    SE_i(a) = spect_exp(a).result.slope(m, t, c, p);
                end
                
                % append
                statement = ['GABA_results(p).rsEEG(m).SE.' time{t} '.' condition{c} ' = SE_i;'];
                eval(statement)      
            end
        end
        
        for c = 1:length(condition)
            % choose data --> SE change
            for a = 1:numel(GABA_results(p).rsEEG(m).SE.bands)
                SE_change_i(a) = spect_exp(a).result.slopechange(m, c, p);
            end

            % append
            statement = ['GABA_results(p).rsEEG(m).SE.change.' condition{c} ' = SE_change_i;'];
            eval(statement)      
        end        
    end
end
clear spect_exp SE_i SE_change_i data_pre statement
clear a c f m p t

% low beta - sigma peak
load('beta.mat'); load('beta_change.mat');
for p = 1:length(participant)
    for m = 1:length(medication)
        GABA_results(p).rsEEG(m).beta.bands = {'beta 1' 'beta 2' 'broad beta'};
        for t = 1:length(time)
            % beta
            statement = ['GABA_results(p).rsEEG(m).beta.' time{t} '(1, 1:3) = squeeze(beta(m, t, p, :));'];
            eval(statement)            
        end
        
        % beta change = normalized
        GABA_results(p).rsEEG(m).beta.change(1, 1:3) = squeeze(beta_change(m, p, :));     
    end
end
clear beta beta_change p m

% save
save('GABA_results.mat', 'GABA_results')

%% 4) extract MEP measures
% load data
load('GABA_MEP.mat')

% extract info 
for p = 1:length(participant)
    for m = 1:length(medication)
        GABA_results(p).MEP(m).medication = medication{m};        
        
        % loop through the output table
        for s = [2,3] % stimulus --> TS, ppTMS
            for t = 1:length(time)
                % choose the row
                rows = (output.subject == p & ...
                    categorical(output.medication) == medication{m} & ...
                    categorical(output.stimulus) == stimulus{s} & ...
                    categorical(output.time) == time{t});
                
                % MEP - with zero reposnses
                statement = ['GABA_results(p).MEP(m).' stimulus{s} '.' time{t} '.amplitude = output.amplitude_zero(rows);'];
                eval(statement)
                statement = ['GABA_results(p).MEP(m).' stimulus{s} '.' time{t} '.epochs = output.epochs_zero(rows);'];
                eval(statement)
                
                % MEP - without zero reposnses
                statement = ['GABA_results(p).MEP(m).' stimulus{s} '.' time{t} '.amplitude_nozero = output.amplitude_nozero(rows);'];
                eval(statement)
                statement = ['GABA_results(p).MEP(m).' stimulus{s} '.' time{t} '.epochs_nozero = output.epochs_nozero(rows);'];
                eval(statement)                               
            end
            
            % calculate MEP change --> % baseline
            rows_pre = (output.subject == p & ...
                    categorical(output.medication) == medication{m} & ...
                    categorical(output.stimulus) == stimulus{s} & ...
                    categorical(output.time) == 'pre');
            rows_post = (output.subject == p & ...
                    categorical(output.medication) == medication{m} & ...
                    categorical(output.stimulus) == stimulus{s} & ...
                    categorical(output.time) == 'post');
            MEP_change = output.amplitude_zero(rows_post) / output.amplitude_zero(rows_pre) * 100;
            statement = ['GABA_results(p).MEP(m).' stimulus{s} '.change = ' num2str(MEP_change) ';'];
            eval(statement)            
        end
        
        % calculate SICI --> % TS
        for t = 1:length(time)
            rows_TS = (output.subject == p & ...
                    categorical(output.medication) == medication{m} & ...
                    categorical(output.stimulus) == 'TS' & ...
                    categorical(output.time) == time{t});
            rows_ppTMS = (output.subject == p & ...
                    categorical(output.medication) == medication{m} & ...
                    categorical(output.stimulus) == 'ppTMS' & ...
                    categorical(output.time) == time{t});
            MEP_SICI = output.amplitude_zero(rows_ppTMS) / output.amplitude_zero(rows_TS) * 100;
            statement = ['GABA_results(p).MEP(m).SICI.' time{t} ' = ' num2str(MEP_SICI) ';'];
            eval(statement)
        end
    end
end
clear output rows statement MEP_change rows_pre rows_post rows_TS rows_ppTMS MEP_SICI
clear p m s t 

% save
save('GABA_results.mat', 'GABA_results')

%% 5) calculate mean values
for m = 1:length(medication)
    for t = 1:length(time)
        if t == 1
            timepoint = 'pre';
        else
            timepoint = 'post';
        end
            
        % rMT
        for p = 1:length(participant)
            statement = ['data(p) = GABA_results(p).session(m).rmt.' timepoint ';'];
            eval(statement)
        end        
        rMT(m, t, 1) = mean(data);
        rMT(m, t, 2) = std(data);
        
        % TEP epochs
        for s = 1:length(stimulus)
            for p = 1:length(participant)                
                statement = ['data(p) = GABA_results(p).TEP(m).epochs.' timepoint '(s);'];
                eval(statement)
            end 
            TEP_epochs(m, t, s, 1) = mean(data);
            TEP_epochs(m, t, s, 2) = std(data);
        end 
        
        % TEP ICA
        for p = 1:length(participant)                
            statement = ['data(p) = GABA_results(p).TEP(m).ICA.' timepoint ';'];
            eval(statement)
        end 
        TEP_ICA(m, t, 1) = mean(data);
        TEP_ICA(m, t, 2) = std(data);
        

        
        % MEP epochs
    end
end
clear m t p s timepoint statement data

% save average values 
save('GABA_average.mat', 'rMT', 'TEP_ICA', 'TEP_epochs')

%% 6) export data for R
% pharmacological activation of GABAARs
% prepare empty table
GABA_med_long = table;
GABA_med_long.subject = zeros(0);  GABA_med_long.medication = zeros(0); GABA_med_long.time = zeros(0); GABA_med_long.stimulus = zeros(0); 
GABA_med_long.peak = zeros(0); GABA_med_long.TEP = zeros(0); GABA_med_long.MEP = zeros(0); 
GABA_med_long.beta = zeros(0); GABA_med_long.AAC = zeros(0); GABA_med_long.SE = zeros(0); 
GABA_med_long.rMT = zeros(0); 

% cycle through entries (= rows)
row_cnt = 1;
for p = 1:length(participant)
    for m = 1:length(medication)
        for t = 1:length(time)
            for s = [1 2]
                for c = 1:length(peaks)
                    % fill corresponding line:
                    % condition info
                    GABA_med_long.subject(row_cnt) = participant(p);             
                    GABA_med_long.medication(row_cnt) = m - 1;
                    GABA_med_long.time(row_cnt) = t - 1;
                    GABA_med_long.stimulus(row_cnt) = s - 1;
                    GABA_med_long.peak(row_cnt) = c;

                    % outcome variables - TEP
                    statement = ['GABA_med_long.TEP(row_cnt) = GABA_results(p).TEP(m).' stimulus{s} '.' time{t} '(c);'];
                    eval(statement)

                    % outcome variables - MEP
                    switch s
                        case 1
                            GABA_med_long.MEP(row_cnt) = 0;
                        case 2
                            statement = ['GABA_med_long.MEP(row_cnt) = GABA_results(p).MEP(m).TS.' time{t} '.amplitude;'];
                            eval(statement)
                    end

                    % RS-EEG: beta 1 - sigma peak
                    statement = ['GABA_med_long.beta(row_cnt) = GABA_results(p).rsEEG(m).beta.' time{t} '(1);'];
                    eval(statement)
                    
                    % RS-EEG: AAC - broad alpha
                    statement = ['GABA_med_long.AAC(row_cnt) = GABA_results(p).rsEEG(m).AAC.' time{t} '(4);'];
                    eval(statement)
                    
                    % RS-EEG: SE - eyes open
                    statement = ['GABA_med_long.SE(row_cnt) = GABA_results(p).rsEEG(m).SE.' time{t} '.open(1);'];
                    eval(statement)

                    % rMT 
                    statement = ['GABA_med_long.rMT(row_cnt) = GABA_results(p).session(m).rmt.' time{t} ';'];
                    eval(statement)

                    % update row count
                    row_cnt = row_cnt + 1;
                end
            end
        end
    end
end
clear row_cnt p m s c statement 

% save as CSV
writetable(GABA_med_long, 'GABA_med_long.csv', 'Delimiter', ',');

%% 7) effect of medication
timepoint = {'1.5h' '2h' '2.5h'};

% extract data
for m = 1:length(medication)
    for p = 1:length(participant)
        for t = 1:length(timepoint)
            data_med(m, p, t) = GABA_results(p).session(m).vigilance(t);
        end
    end
end
clear m p t

% extract average values
for m = 1:length(medication)
    for t = 1:length(timepoint)
        data_med_avg(m, t, 1) = mean(data_med(m, :, t), 'omitnan');
        data_med_avg(m, t, 2) = std(data_med(m, :, t), 'omitnan');          % STD
        n = sum(~isnan(data_med(m, :, t)));
        data_med_avg(m, t, 3) = data_med_avg(m, t, 2)/sqrt(n);              % SEM
        data_med_avg(m, t, 4) = data_med_avg(m, t, 3) * z;                  % 95% confidence interval      
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
    err(m) = errorbar(x, data_med_avg(m, :, 1), data_med_avg(m, :, 4));

    % add parameters
    err(m).Color = colours2((m - 1)*2 + 2, :);
    err(m).LineWidth = 2;
    err(m).Marker = 'o';
    err(m).MarkerFaceColor = colours2((m - 1)*2 + 2, :);
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
savefig([fig_name '.fig'])
saveas(fig, [fig_name '.png'])

% update counter
figure_counter = figure_counter + 1;

% ----- test for differences - rm-ANOVA -----
% prepare data table
data_table = table; data_table.participant = participant'; 
for m = 1:length(medication)
    for t = 1:length(timepoint)
        statement = ['data_table.' medication{m} '_t' num2str(t) ' = data_med(m, :, t)'';'];
        eval(statement)
    end
end
clear m t statement 
 
% within-subjects design
wd = table([repmat(medication(1), length(timepoint), 1); repmat(medication(2), length(timepoint), 1)],...
    [1:length(timepoint) 1:length(timepoint)]','VariableNames',{'treatment','timepoint'});
wd.treatment = categorical(wd.treatment); wd.timepoint = categorical(wd.timepoint);

% rm-ANOVA
rm = fitrm(data_table, 'placebo_t1-alprazolam_t3~1', 'WithinDesign', wd);

% check for approx. normal distribution
figure(figure_counter)
histogram(reshape(data_table{:,2:end},[],1), 10)

% check for sphericity
rm.mauchly

% results table
ranova(rm, 'WithinModel', 'treatment*timepoint') 
figure(figure_counter)
boxplot(data_table{:,2:end})
xlabel('measurements')
ylabel('arousal assessment')

%% ) TEP x RS-EEG correlation
% parameters
varnames = [peaks, {'AAC' 'SEo' 'SEc' 'beta'}];

% ----- extract data -----
for m = 1:length(medication)
    for s = [1 2]
        % change in TEPs 
        for k = 1:length(peaks)
            for p = 1:numel(GABA_results)                
                statement = ['data_cor(m, s, p, k) = GABA_results(p).TEP(m).' stimulus{s} '.change(k);'];            
                eval(statement)    
            end
        end
        
        % change in AAC - high alpha 3 only
        for p = 1:numel(GABA_results)            
            data_cor(m, s, p, k+1) = GABA_results(p).rsEEG(m).AAC.change(3);              
        end
        
        % change in SE - lower frequencies only
        for p = 1:numel(GABA_results)            
            data_cor(m, s, p, k+2) = GABA_results(p).rsEEG(m).SE.change.open(1);   
            data_cor(m, s, p, k+3) = GABA_results(p).rsEEG(m).SE.change.closed(1);  
        end
        
        % broad beta 
        for p = 1:numel(GABA_results)            
            data_cor(m, s, p, k+4) = GABA_results(p).rsEEG(m).beta.change;     
        end
    end
end
clear statement m s p k

% ----- plot data to assess linear correlation between all variables -----
for m = 1:length(medication)
    for s = [1 2]
        % select the data                
        mat_cor = squeeze(data_cor(m, s, :, :));
        
        % create correlation matrix
        fig = correlation_preview(mat_cor, varnames, 'method', 'Pearson'); 
        
        % save figure
        figure_name = ['corr_NORMAL_overview_' medication{m} '_' stimulus{s}];
        savefig([figure_name])
        saveas(fig, [figure_name '.png'])             
    end
end
clear m s fig figure_name

% ----- significant linear correlation: sessions separately -----
for m = 1:length(medication)
    for s = [1 2]          
        % select the data                
        mat_cor = squeeze(data_cor(m, s, :, :));
        
        % identify significant cases        
        [cor_coef, cor_p] = corrcoef(mat_cor);
        [row, col] = find(cor_p < alpha);
        
        % plot significant cases
        for a = 1:length(row)    
            % prepare linear model: y ~ 1 + x
            data_model = fitlm(mat_cor(:, row(a)), mat_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);
            
            % choose only correlations that show TEP-MEP interactions
            if ismember(col(a), 6:10)            
                % plot data + regression line
                fig = figure(figure_counter);
                hold on
                plot_cor = plotAdded(data_model);

                % adjust parameters    
                title([medication{m} ' - ' stimulus{s} ' : ' varnames{col(a)} ' ~ ' varnames{row(a)}])
                xlabel(['change in ' varnames{row(a)}]); ylabel(['change in ' varnames{col(a)}]);
                set(gca, 'FontSize', 14)
                plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
                plot_cor(1).MarkerEdgeColor = colours2(2, :); plot_cor(1).MarkerFaceColor = colours2(2, :);
                plot_cor(2).Color = colours2(3, :); plot_cor(2).LineWidth = 2; 
                plot_cor(3).Color = colours2(3, :); plot_cor(3).LineWidth = 2;
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
                set(T(3), 'fontsize', 14, 'color', colours2(3, :)); 

                % save figure and continue
                figure_name = ['corr_NORMAL_' medication{m} '_' stimulus{s} '_' varnames{row(a)} '_' varnames{col(a)}];
                savefig([figure_name])
                saveas(fig, [figure_name '.png'])
                figure_counter = figure_counter + 1;
            end
        end
    end
end
clear mat_cor cor_coef cor_p data_model plot_cor row col fig figure_name T text_pos m s a  

% ----- plot data to assess non-linear correlation between all variables -----
for m = 1:length(medication)
    for s = [1 2]
        % select the data                
        mat_cor = squeeze(data_cor(m, s, :, :));
        
        % create correlation matrix
        fig = correlation_preview(mat_cor, varnames, 'method', 'Spearman'); 
        
        % save figure
        figure_name = ['corr_RANKED_overview_' medication{m} '_' stimulus{s}];
        savefig([figure_name])
        saveas(fig, [figure_name '.png'])             
    end
end
clear m s fig figure_name

% ----- significant non-linear correlation: sessions separately -----
for m = 1:length(medication)
    for s = 2          
        % select the data                
        mat_cor = squeeze(data_cor(m, s, :, :));
        
        % identify significant cases        
        [cor_coef, cor_p] = corr(mat_cor, 'Type', 'Spearman');
        [row, col] = find(cor_p < alpha);
        
        % rank the data
        for a = 1:size(mat_cor, 2)
            [temp, mat_cor(:, a)]  = ismember(mat_cor(:, a), unique(mat_cor(:, a)));
        end
        
        % plot significant cases
        for a = 1:length(row)    
            % prepare linear model: y ~ 1 + x
            data_model = fitlm(mat_cor(:, row(a)), mat_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);
            
            % choose only correlations that show TEP-MEP interactions
            if ismember(col(a), 6:10)            
                % plot data + regression line
                fig = figure(figure_counter);
                hold on
                plot_cor = plotAdded(data_model);

                % adjust parameters    
                title([medication{m} ' - ' stimulus{s} ' : ' varnames{col(a)} ' ~ ' varnames{row(a)}])
                xlabel(['change in ' varnames{row(a)}]); ylabel(['change in ' varnames{col(a)}]);
                set(gca, 'FontSize', 14)
                plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
                plot_cor(1).MarkerEdgeColor = colours2(2, :); plot_cor(1).MarkerFaceColor = colours2(2, :);
                plot_cor(2).Color = colours2(3, :); plot_cor(2).LineWidth = 2; 
                plot_cor(3).Color = colours2(3, :); plot_cor(3).LineWidth = 2;
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
                set(T(3), 'fontsize', 14, 'color', colours2(3, :)); 

                % save figure and continue
                figure_name = ['corr_RANKED_' medication{m} '_' stimulus{s} '_' varnames{row(a)} '_' varnames{col(a)}];
                savefig([figure_name])
                saveas(fig, [figure_name '.png'])
                figure_counter = figure_counter + 1;
            end
        end
    end
end
clear varnames data_cor mat_cor cor_coef cor_p data_model plot_cor row col fig figure_name T text_pos m s a  


%% ) TEP x MEP correlation
% parameters
varnames = [peaks, {'MEP'}];

% ----- extract data -----
for m = 1:length(medication)
    % change in TEPs 
    for k = 1:length(peaks)
        for p = 1:numel(GABA_results)                
            data_cor(m, p, k) = GABA_results(p).TEP(m).TS.change(k);             
        end
    end

    % change in MEPs
    for p = 1:numel(GABA_results)            
        data_cor(m, p, k+1) = GABA_results(p).MEP(m).TS.change;              
    end
end
clear m p k

% ----- significant linear correlation: sessions separately -----
for m = 1:length(medication)
    % calculate correlation coefficient and p values
    mat_cor = squeeze(data_cor(m, :, :));
    [cor_coef, cor_p] = corrcoef(mat_cor);

    % identify significant cases 
    [row, col] = find(cor_p < alpha );

    % plot significant cases
    for a = 1:length(row)    
        % prepare linear model: y ~ 1 + x
        data_model = fitlm(mat_cor(:, row(a)), mat_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);
        
        % choose only correlations that show TEP-MEP interactions
        if col(a) == 6
            % plot data + regression line
            fig = figure(figure_counter);
            hold on
            plot_cor = plotAdded(data_model);

            % adjust parameters    
            title([medication{m} ' : ' varnames{col(a)} ' ~ ' varnames{row(a)}])
            xlabel(['change in ' varnames{row(a)}]); ylabel(['change in ' varnames{col(a)}]);
            set(gca, 'FontSize', 14)
            plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
            plot_cor(1).MarkerEdgeColor = colours2(2, :); plot_cor(1).MarkerFaceColor = colours2(2, :);
            plot_cor(2).Color = colours2(3, :); plot_cor(2).LineWidth = 2; 
            plot_cor(3).Color = colours2(3, :); plot_cor(3).LineWidth = 2;
            legend off
            if data_model.Coefficients.Estimate(2) > 0
                text_pos = [0.95 0.85 0.75];
            else
                text_pos = [0.25 0.15 0.05];
            end
            T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
            T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
            T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.3f', cor_coef(row(a), col(a)), cor_p(row(a), col(a))), 'Units', 'Normalized');
            set(T(1), 'fontsize', 14, 'fontweight', 'bold', 'fontangle', 'italic');   
            set(T(2), 'fontsize', 14); 
            set(T(3), 'fontsize', 14, 'color', colours2(3, :)); 
            hold off

            % save figure and continue
            figure_name = ['corr_NORMAL_' medication{m} '_' varnames{row(a)} '_' varnames{col(a)}];
            savefig([figure_name])
            saveas(fig, [figure_name '.png'])
            figure_counter = figure_counter + 1;
        end
    end
end
clear mat_cor cor_coef cor_p data_model plot_cor row col fig figure_name T text_pos m a    

% ----- significant non-linear correlation: sessions separately -----
% plot significant cases
for m = 1:length(medication)
    % calculate correlation coefficient and p values
    mat_cor = squeeze(data_cor(m, :, :));
    [cor_coef, cor_p] = corr(mat_cor, 'Type', 'Spearman');
    
    % identify significant cases 
    [row, col] = find(cor_p < alpha);
    
    % rank the data
    for a = 1:size(mat_cor, 2)
        [temp, mat_ranked(:, a)]  = ismember(squeeze(mat_cor(:, a)), unique(squeeze(mat_cor(:, a))));
    end
    
    % plot correlations
    for a = 1:length(row)        
        % prepare linear model: y ~ 1 + x        
        data_model = fitlm(mat_ranked(:, row(a)), mat_ranked(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

        % choose only correlations that show TEP-MEP interactions
        if col(a) == 6
            % plot data + regression line
            fig = figure(figure_counter);
            hold on
            plot_cor = plotAdded(data_model);

            % adjust parameters    
            title([medication{m} ' : ' varnames{col(a)} ' ~ ' varnames{row(a)}])
            xlabel(['change in ' varnames{row(a)} ' - ranked']); ylabel(['change in ' varnames{col(a)} ' - ranked']);
            set(gca, 'FontSize', 14)
            plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
            plot_cor(1).MarkerEdgeColor = colours2(2, :); plot_cor(1).MarkerFaceColor = colours2(2, :);
            plot_cor(2).Color = colours2(3, :); plot_cor(2).LineWidth = 2; 
            plot_cor(3).Color = colours2(3, :); plot_cor(3).LineWidth = 2;
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
            set(T(3), 'fontsize', 14, 'color', colours2(3, :)); 
            hold off

            % save figure and continue
            figure_name = ['corr_RANKED_' medication{m} '_' varnames{row(a)} '_' varnames{col(a)}];
            savefig([figure_name])
            saveas(fig, [figure_name '.png'])
            figure_counter = figure_counter + 1;
        end
    end
end
clear varnames data_cor mat_cor cor_coef cor_p data_model plot_cor row col fig figure_name T text_pos a   

%% ) TEP SICI x MEP SICI correlation
% parameters
varnames = [peaks, {'MEP'}];

% ----- extract data -----
for m = 1:length(medication)
    % change in TEPs 
    for k = 1:length(peaks)
        for p = 1:numel(GABA_results)                
            data_cor(m, p, k) = GABA_results(p).TEP(m).SICI.pre(k);             
        end
    end

    % change in MEPs
    for p = 1:numel(GABA_results)            
        data_cor(m, p, k+1) = GABA_results(p).MEP(m).SICI.pre;              
    end
end
clear m p k

% ----- significant linear correlation: sessions separately -----
% plot significant cases
for m = 1:length(medication)
    % calculate correlation coefficient and p values
    mat_cor = squeeze(data_cor(m, :, :));
    [cor_coef, cor_p] = corrcoef(mat_cor);
    
    % identify significant cases 
    [row, col] = find(cor_p < alpha);
    
    % plot correlations
    for a = 1:length(row) 
        % prepare linear model: y ~ 1 + x        
        data_model = fitlm(mat_cor(:, row(a)), mat_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

        % choose only correlations that show TEP-MEP interactions
        if col(a) == 6
            % plot data + regression line
            fig = figure(figure_counter);
            hold on
            plot_cor = plotAdded(data_model);

            % adjust parameters    
            title([medication{m} ' - SICI : ' varnames{col(a)} ' ~ ' varnames{row(a)}])
            xlabel(['SICI in ' varnames{row(a)}]); ylabel(['SICI in ' varnames{col(a)}]);
            set(gca, 'FontSize', 14)
            plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
            plot_cor(1).MarkerEdgeColor = colours2(2, :); plot_cor(1).MarkerFaceColor = colours2(2, :);
            plot_cor(2).Color = colours2(3, :); plot_cor(2).LineWidth = 2; 
            plot_cor(3).Color = colours2(3, :); plot_cor(3).LineWidth = 2;
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
            set(T(3), 'fontsize', 14, 'color', colours2(3, :)); 
            hold off

            % save figure and continue
            figure_name = ['corr_NORMAL_' medication{m} '_SICI_' varnames{row(a)} '_' varnames{col(a)}];
            savefig([figure_name])
            saveas(fig, [figure_name '.png'])
            figure_counter = figure_counter + 1;
        end
    end
end
clear mat_cor cor_coef cor_p data_model plot_cor row col fig figure_name T text_pos a   

% ----- significant non-linear correlation: sessions separately -----
% plot significant cases
for m = 1:length(medication)
    % calculate correlation coefficient and p values
    mat_cor = squeeze(data_cor(m, :, :));
    [cor_coef, cor_p] = corr(mat_cor, 'Type', 'Spearman');
    
    % identify significant cases 
    [row, col] = find(cor_p < alpha);
    
    % rank the data
    for a = 1:size(mat_cor, 2)
        [temp, mat_ranked(:, a)]  = ismember(squeeze(mat_cor(:, a)), unique(squeeze(mat_cor(:, a))));
    end
    
    % plot correlations
    for a = 1:length(row)        
        % prepare linear model: y ~ 1 + x        
        data_model = fitlm(mat_ranked(:, row(a)), mat_ranked(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

        % choose only correlations that show TEP-MEP interactions
        if col(a) == 6
            % plot data + regression line
            fig = figure(figure_counter);
            hold on
            plot_cor = plotAdded(data_model);

            % adjust parameters    
            title([medication{m} ' - SICI : ' varnames{col(a)} ' ~ ' varnames{row(a)}])
            xlabel(['SICI in ' varnames{row(a)} ' - ranked']); ylabel(['SICI in ' varnames{col(a)} ' - ranked']);
            set(gca, 'FontSize', 14)
            plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
            plot_cor(1).MarkerEdgeColor = colours2(2, :); plot_cor(1).MarkerFaceColor = colours2(2, :);
            plot_cor(2).Color = colours2(3, :); plot_cor(2).LineWidth = 2; 
            plot_cor(3).Color = colours2(3, :); plot_cor(3).LineWidth = 2;
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
            set(T(3), 'fontsize', 14, 'color', colours2(3, :)); 
            hold off

            % save figure and continue
            figure_name = ['corr_RANKED_' medication{m} '_SICI_' varnames{row(a)} '_' varnames{col(a)}];
            savefig([figure_name])
            saveas(fig, [figure_name '.png'])
            figure_counter = figure_counter + 1;
        end
    end
end
clear varnames data_cor mat_cor cor_coef cor_p data_model plot_cor row col fig figure_name T text_pos a   

%% ) TEP change x TEP SICI correlation
% parameters
for a = 1:length(peaks)
    peaks_SICI{a} = [peaks{a}, ' SICI'];
end
varnames = [peaks, peaks_SICI];

 % ----- extract data -----
for m = 1:length(medication)
    % change in TEPs 
    for k = 1:length(peaks)
        for p = 1:numel(GABA_results)                
            statement = ['data_cor(m, p, k) = GABA_results(p).TEP(m).TS.change(k);'];            
            eval(statement)    
        end
    end

    % TEP SICI
    for k = 1:length(peaks)
        for p = 1:numel(GABA_results)                
            data_cor(m, p, length(peaks) + k) = GABA_results(p).TEP(m).SICI.pre(k);              
        end
    end
end
clear peaks_SICI a m p k

% ----- plot data to assess linear correlation between all variables -----
for m = 1:length(medication)
    % select the data                
    mat_cor = squeeze(data_cor(m, :, :));

    % create correlation matrix
    fig = correlation_preview(mat_cor, varnames, 'method', 'Pearson'); 

    % save figure
    figure_name = ['corr_NORMAL_overview_' medication{m} '_TEP-SICI'];
    savefig([figure_name])
    saveas(fig, [figure_name '.png'])             
end
clear m s fig figure_name

% ----- significant linear correlation: only alprazolam -----
% calculate correlation coefficient and p values
mat_cor = squeeze(data_cor(2, :, :));
[cor_coef, cor_p] = corr(mat_cor, 'Type', 'Pearson');

% identify significant cases 
[row, col] = find(cor_p < alpha);

% plot significant cases
for a = 1:length(row)  
    if (row(a) <= length(peaks) & col(a) > length(peaks)) | (row(a) > length(peaks) & col(a) <= length(peaks))
        if (row(a) == length(peaks) + col(a)) | (col(a) == length(peaks) + row(a))
            % prepare linear model: y ~ 1 + x
            data_model = fitlm(mat_cor(:, row(a)), mat_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

            % plot data + regression line
            fig = figure(figure_counter);
            hold on
            plot_cor = plotAdded(data_model);

            % adjust parameters    
            title(['alprazolam : ' varnames{row(a)} ' ~ ' varnames{col(a)}])
            xlabel(['change in ' varnames{row(a)}]); ylabel(['change in ' varnames{col(a)}]);
            set(gca, 'FontSize', 14)
            plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
            plot_cor(1).MarkerEdgeColor = colours2(2, :); plot_cor(1).MarkerFaceColor = colours2(2, :);
            plot_cor(2).Color = colours2(3, :); plot_cor(2).LineWidth = 2; 
            plot_cor(3).Color = colours2(3, :); plot_cor(3).LineWidth = 2;
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
            set(T(3), 'fontsize', 14, 'color', colours2(3, :));   
            hold off

            % save figure and continue
            figure_name = ['corr_NORMAL_alprazolam_TEP-SICI_' varnames{row(a)} '_' varnames{col(a)}];
            savefig([figure_name])
            saveas(fig, [figure_name '.png'])
            figure_counter = figure_counter + 1;
        end
    end
end
clear mat_cor cor_coef cor_p data_model plot_cor row col fig figure_name T text_pos m a 

% ----- plot data to assess non-linear correlation between all variables -----
for m = 1:length(medication)
    % select the data                
    mat_cor = squeeze(data_cor(m, :, :));

    % create correlation matrix
    fig = correlation_preview(mat_cor, varnames, 'method', 'Spearman'); 

    % save figure
    figure_name = ['corr_RANKED_overview_' medication{m} '_TEP-SICI'];
    savefig([figure_name])
    saveas(fig, [figure_name '.png'])             
end
clear m s fig figure_name

% ----- significant linear correlation: only alprazolam -----
% calculate correlation coefficient and p values
mat_cor = squeeze(data_cor(2, :, :));
[cor_coef, cor_p] = corr(mat_cor, 'Type', 'Spearman');

% identify significant cases 
[row, col] = find(cor_p < alpha);

% rank the data
for a = 1:size(mat_cor, 2)
    [temp, mat_cor(:, a)]  = ismember(mat_cor(:, a), unique(mat_cor(:, a)));
end

% plot significant cases
for a = 1:length(row)  
    if (row(a) <= length(peaks) & col(a) > length(peaks)) | (row(a) > length(peaks) & col(a) <= length(peaks))
        if (row(a) == length(peaks) + col(a)) | (col(a) == length(peaks) + row(a))
            % prepare linear model: y ~ 1 + x
            data_model = fitlm(mat_cor(:, row(a)), mat_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

            % plot data + regression line
            fig = figure(figure_counter);
            hold on
            plot_cor = plotAdded(data_model);

            % adjust parameters    
            title(['alprazolam : ' varnames{row(a)} ' ~ ' varnames{col(a)}])
            xlabel(['change in ' varnames{row(a)}]); ylabel(['change in ' varnames{col(a)}]);
            set(gca, 'FontSize', 14)
            plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
            plot_cor(1).MarkerEdgeColor = colours2(2, :); plot_cor(1).MarkerFaceColor = colours2(2, :);
            plot_cor(2).Color = colours2(3, :); plot_cor(2).LineWidth = 2; 
            plot_cor(3).Color = colours2(3, :); plot_cor(3).LineWidth = 2;
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
            set(T(3), 'fontsize', 14, 'color', colours2(3, :));   
            hold off

            % save figure and continue
            figure_name = ['corr_RANKED_alprazolam_TEP-SICI_' varnames{row(a)} '_' varnames{col(a)}];
            savefig([figure_name])
            saveas(fig, [figure_name '.png'])
            figure_counter = figure_counter + 1;
        end
    end
end
clear varnames data_cor mat_cor cor_coef cor_p data_model plot_cor row col fig figure_name T text_pos m a 



