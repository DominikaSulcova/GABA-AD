function  filename = GABA_initialize_logfile(session_info)
% ------------------------------------------------------------------------
% Fnc: Creates a text logfile that serves as documentation for GABA-AD project 
% Input: session_info (cell array)
% Author: Dominika (2020) 
% ------------------------------------------------------------------------

% define the filename
filename = ['GABA-AD_' char(session_info{1}) char(session_info{2}) '.txt'];

% write the file 
fileID = fopen(filename,'w');
fprintf(fileID, '******************************************************************************************************\r\n');
fprintf(fileID, 'study: GABA-AD\r\n'); 
fprintf(fileID, ['subject: ' [char(session_info{1}) char(session_info{2})] '\r\n']);
fprintf(fileID, ['recorded: ' char(session_info{7}) '\r\n']);
fprintf(fileID, 'author: Dominika Sulcova\r\n');
fprintf(fileID, '******************************************************************************************************\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, 'DATA ACQUISITION\r\n');
fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, '- two recording blocks - pre and post medication (alprazolam  1mg orally)\r\n');
fprintf(fileID, '       saccadic peak velocity measurement (Pupil wearable camera system)\r\n');
fprintf(fileID, '       RS-EEG 3 mins eyes open + 3 mins eyes closed\r\n');
fprintf(fileID, '       TMS-EEG-EMG over M1\r\n');
fprintf(fileID, '           3 conditions per 60 stim - spTMS_TS(120 %%rMT); spTMS_CS(80 %%rMT); ppTMS(80+120 %%rMT; ISI 2.5ms)\r\n');
fprintf(fileID, '           --> 3 blocks of 60 stim, coditions randomly mixed\r\n');
fprintf(fileID, '       TMS-EEG-EMG over AG\r\n');
fprintf(fileID, '           80 stim, 120 %%rMT\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '- TMS:     intensity in %%MSO, rMT has was not adjusted in the post-medication block\r\n');
fprintf(fileID, '           MagVenture MagPro X100 (biphasic sin pulse) + Visor2 neuronavigation (generic brain model)\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '- EEG:     Bittium system - TESLA amplifiers, NeurOne recording software, 20kHz sampling rate, 3500Hz LP device filter\r\n');
fprintf(fileID, '           setup with 32 electrodes (10-20), from which M1(29) + M2(30) were used as references, ground at the AFz position\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '- EMG:     recorded from contralateral FDI muscle\r\n');
fprintf(fileID, '           MOBI system connected to Visor2, 1024Hz sampling rate\r\n');
fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, 'beginning:\r\n'); 
fprintf(fileID, 'medication:\r\n'); 
fprintf(fileID, 'end:\r\n'); 
fprintf(fileID, '\r\n');
fprintf(fileID, 'rMT pre:\r\n'); 
fprintf(fileID, 'rMT post:\r\n'); 
fprintf(fileID, '\r\n');
fprintf(fileID, 'subjective 1.5h:\r\n'); 
fprintf(fileID, 'subjective 2.0h:\r\n'); 
fprintf(fileID, 'subjective 2.5h:\r\n'); 
fprintf(fileID, '\r\n');
fprintf(fileID, 'closest electrodes M1:\r\n'); 
fprintf(fileID, 'PRE - \r\n'); 
fprintf(fileID, 'POST - \r\n'); 
fprintf(fileID, 'closest electrodes AG:\r\n'); 
fprintf(fileID, 'PRE - \r\n'); 
fprintf(fileID, 'POST - \r\n'); 
fprintf(fileID, '\r\n');
fprintf(fileID, 'notes:\r\n'); 
fprintf(fileID, 'M1 blockX/AG - EMG xxx - OUT \r\n'); 
fprintf(fileID, '\r\n');
fprintf(fileID, '\r\n');
fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
fprintf(fileID, 'DATA IMPORT\r\n');
fprintf(fileID, '------------------------------------------------------------------------------------------------------\r\n');
fclose(fileID);
end