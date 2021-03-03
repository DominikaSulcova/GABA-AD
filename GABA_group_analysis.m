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
%       - mean amplitude over fbands
%       - amplitude change
%       - alpha attenuation coeficient
%       - spectral exponent beta 
% 4) MEP
% 
% ----- STATISTICS -----
% 5)

%% parameters
% dataset
participant = 1:20;
stimulus = {'CS' 'TS' 'ppTMS'};
medication = {'placebo' 'alprazolam'};
time = {'pre' 'post'};

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
        for a = 1:3 
            if strcmp(class(params{p, 14 + a}), 'cell')
                GABA_results(p).session(1).vigilance(a) = 0;                    % 0 = missing value
            else
                GABA_results(p).session(1).vigilance(a) = params{p, 14 + a};
            end
        end

        GABA_results(p).session(2).medication = 'alprazolam';
        GABA_results(p).session(2).rmt.pre = params{p, 'xanax_rMT_pre'}; 
        GABA_results(p).session(2).rmt.post = params{p, 'xanax_rMT_post'}; 
        GABA_results(p).session(2).rmt.change = params{p, 'xanax_rMT_post'} - params{p, 'xanax_rMT_pre'}; 
        for a = 1:3 
            if strcmp(class(params{p, 7 + a}), 'cell')
                GABA_results(p).session(2).vigilance(a) = 0;
            else
                GABA_results(p).session(2).vigilance(a) = params{p, 7 + a};
            end
        end
    end
else
    load('GABA_results.mat')
end

%% 2) extract TEP measures
% fill in TEP related info
for p = 1:length(participant)
    GABA_results(p).TEP(1).medication = 'placebo';
    GABA_results(p).TEP(1).ICA.pre = params{p, 'zyrtec_ICA_pre'};
    GABA_results(p).TEP(1).ICA.post = params{p, 'zyrtec_ICA_post'};
    
    GABA_results(p).TEP(2).medication = 'alprazolam';
    GABA_results(p).TEP(2).ICA.pre = params{p, 'xanax_ICA_pre'};
    GABA_results(p).TEP(2).ICA.post = params{p, 'xanax_ICA_post'};
end
clear params

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
% dataset info --> ICs, epochs etc.
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

% spectral exponent beta 
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

% save
save('GABA_results.mat', 'GABA_results')

%% 4) extract MEP measures

%% 5) TEP/RS-EEG correlation
% parameters
peaks = {'P30' 'N45' 'P60' 'N100' 'P180'};
varnames = [peaks, {'AAC' 'SEo' 'SEc'}];

% extract data
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
    end
end
clear statement m s p k

% plot correlation matrix
varnames = [peaks, {'AAC' 'SEo' 'SEc'}];
for m = 1:length(medication)
    for s = [1 2]        
        corrplot(squeeze(data_cor(m, s, :, :)), 'type', 'Spearman', 'testR', 'on', 'varNames', varnames);
        set(gcf, 'Name', [medication{m} ' - ' stimulus{s} ' TEP change x RS-EEG marker'])
        savefig(gcf, ['correlation_' medication{m} '_' stimulus{s} '_RS-EEG.fig'])
    end
end
clear varnames data_cor m s     

%% 5) SICI/RS-EEG correlation
% extract data
for m = 1:length(medication)
    % change in SICI
    for k = 1:length(peaks)
        for p = 1:numel(GABA_results)                
            data_cor(m, p, k) = GABA_results(p).TEP(m).SICI.pre(k);               
        end
    end

    % change in AAC - high alpha 3 only
    for p = 1:numel(GABA_results)            
        data_cor(m, p, k+1) = GABA_results(p).rsEEG(m).AAC.change(3);              
    end

    % change in SE - lower frequencies only
    for p = 1:numel(GABA_results)            
        data_cor(m, p, k+2) = GABA_results(p).rsEEG(m).SE.change.open(1);   
        data_cor(m, p, k+3) = GABA_results(p).rsEEG(m).SE.change.closed(1);  
    end
end
clear m p k

% plot correlation matrix
for m = 1:length(medication)    
    corrplot(squeeze(data_cor(m, :, :)), 'type', 'Spearman', 'testR', 'on', 'varNames', varnames);
    set(gcf, 'Name', [medication{m} ' - SICI x RS-EEG markers'])
    savefig(gcf, ['correlation_' medication{m} '_SICI_RS-EEG.fig'])
end
clear data_cor m 

















