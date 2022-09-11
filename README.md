# GABA-AD
PhD project

------------------------------------------------------------
PROJECT DESCRIPTION
------------------------------------------------------------

Experimental questions
------------------------------------------------------------
The aims of the project are:
1)  to characterize the inhibitory neurotransmission mediated by GABA A receptors in young healthy human subjects
    using TMS-EEG/EMG and resting state EEG (RS-EEG)
2)  to identify EEG/EMG biomarkers most suited to assess the state of GABAergic neurotransmission for future comparison of
    patient groups

Two different types of GABAA-mediated inhibition were investigated:
1)  global activation of the GABAergic inhibitory system induced by GABA A receptor stimulation with benzodiazepines
    --> changes observed in TEP and MEP amplitude and in RS-EEG variables
2)  local, phasic GABA A activation induced by sub-threshold conditioning TMS pulse (CS) preceding the testing stimulus (TS)
    --> short-latency intracortical inhibition (SICI) = amplitude reduction of TEPs and MEPs

_**Question 1: How does pharmacological stimulation of GABAARs with alprazolam influence the outcome measures?**
--> Is there a change in TEPs, MEPs or RS-EEG that can be only observed after alprazolam, and not placebo?
--> Is there a correlation across alprazolam-induced changes in different measures?_

_**Question 2: Can we observe changes in TEP and MEP amplitudes following the paired pulse protocol?**
--> Can we identify TEP component(s) that predict the MEP SICI?
--> How does TEP SICI relate to the previously described effect of alprazolam?_

Experimental protocol
------------------------------------------------------------
Datasets for each subject:
- young controls: 2 drug sessions --> alprazolam (Xanax 1mg) + placebo (zyrtec)
- aged subjects: 1 drug session --> alprazolam (xanax 1mg)
- 2 timepoint per session: before + after drug administration
- 2 blocks of RS-EEG per timepoint: 3 mins eyes open + 3 mins eyes closed
- 4 stimulation blocks per timepoint: 3 blocks TMS over left M1 + 1 block TMS over left AG

TMS:
- MagVenture MagPro X100, biphasic sin pulse
- 60 stimuli per stimulation block
- over M1:
  - 3 types of stimuli: 60x spTMS 120 %rMT (TS) + 60x spTMS 80 %rMT (CS) + 60x ppTMS 80/120 %rMT with 2.5 ms ISI (ppTMS) - SICI
  - randomly delivered with variable ISI (4 - 6 s)
- over AG:
	- young controls: 60x spTMS 100 %rMT
  - aged subjects: 80x spTMS 120 %rMT

EEG:
- NeurOne EEG system (MEGA/Bittium)
- 32 (30) electrodes – referenced to mastoids (ground an extra electrode on the front)
- SR 20 kHz

EMG:
- Visor2 MobiEMG amplifier
- SR 1024 Hz

Data analysis
------------------------------------------------------------
There were several types of data collected within the project and analysed using different approaches. In addition,
resting motor threshold (rMT) and subjective vigilance assessment were included in outcome variables.  

1)  Resting state EEG was recorded in order to perform time-frequency transform using simple FFT and amplitudes
    were calculated across available frequencies (up to 45Hz). The resulting spectrum then served to extract
    4 outcome variables:
    - beta increase
    - alpha attenuation coefficient
        --> alpha amplitude eyes closed /  alpha amplitude eyes open
        --> indicates the level of sleepiness, subjective loss of alertness
    - spectral exponent (separately for eyes open and closed)
        --> describes the slope of a least-squares line that can be fitted to the log-transformed power spectrum density
        --> indicates objective deterioration of consciousness, reflects the shift in excitation/inhibition ratio

2)  TMS-evoked potentials (TEPs) were recorded following TMS stimulation of left primary motor cortex (M1) and a non-motor
    cortical area (left dorsal angular gyrus - AG). Simple TMS pulses of sub- and supra-threshold intensities were applied to M1,
    AG was stimulated with 100%rMT. TEPs were also recorded following the application of paired-pulse protocol over the M1 aiming
    to record SICI.

    Following data preprocessing and cleaning, TEPs were characterised at the group level using **Microstate analysis** as implemented
    in the RAGU software in order to identify the duration and topographic distribution of individual TEP components. Based on these,
    3 electrodes of interest (EOIs) located at the maximum amplitude were chosen for each component, and the pooled signal from these EOIs
    was used to extract peak amplitude for each individual subject and condition.
    Rather than taking one-datapoint peak amplitude, windows of default lengths were first determined for each TEP component
    and their latency was individually adjusted for each dataset, then the mean amplitude was calculated as the average
    of 25% of most prominent datapoints within the window. Peak latency was calculated as the median x value of these averaged datapoints.

3)  Motor-evoked potentials (MEPs) were recorded at the same time as TEPs, clear motor responses were obtained following TS and ppTMS
    stimulation (not CS). Trials containing baseline muscular activity were semi-automatically discarded, trials showing no response
    were left in the dataset. Peak-to-peak MEP amplitude was then calculated for each subject and condition.

Statistical analysis

Detailed, step-by-step descriptions of data preprocessing, extraction of outcome variables and subsequent statistical analysis
are available in the repository as following word documents:
- GABA-AD_PROTOCOL_RS-EEG
- GABA-AD_PROTOCOL_TMS-EEG
- GABA-AD_PROTOCOL_TMS-EMG

------------------------------------------------------------
AVAILABLE SCRIPTS
------------------------------------------------------------

RS-EEG
------------------------------------------------------------
preprocessing:
- GABA_rsEEG_preprocessing1.m
- GABA_rsEEG_preprocessing2.m
- GABA_rsEEG_preprocess.lwscript
- GABA_rsEEG_umica.lwscript

RS-EEG analysis:
- GABA_rsEEG_amplitude.m
- GABA_rsEEG_SE.m
- GABA_rsEEG_YC_process.m
- GABA_rsEEG_visualization.m


EEG: TEPs from M1 and AG
------------------------------------------------------------
import:
- GABA_EEG_import.m
- EEG_import_MEGA.m + EEG_history_import.mat
- GABA_initialize logfile.m

preprocessing:
- GABA_EEG_parse_events.m
- GABA_EEG_unblind.m
- GABA_EEG_maintenance.m
- GABA_EEG_correct_ppTMS.m
- GABA_EEG_preprocess.lwscript
- GABA_EEG_ffilt_crop.lwscript
- GABA_EEG_ica_FFT.lwscript
- GABA_EEG_ica_timecourse.lwscript

TEP analysis:
- GABA_YC_topography.m
- GABA_TEP_YC_M1_process.m
- GABA_TEP_YC_AG_process.m

EMG: MEPs from M1
------------------------------------------------------------
import:
- GABA_EMG_import.m
- EMG_import_VHDR.m + EMG_history_import.mat

preprocessing + MEP analysis:
- GABA_EMG_process.m

Group analysis & statistics
------------------------------------------------------------
- GABA_YC_export.m
- GABA_YC_stats_medication.rmd
- GABA_YC_stats_SICI.rmd
- GABA_YC_stats_SICIxmedication.rmd
