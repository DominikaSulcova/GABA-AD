# GABA-AD
PhD project

Experimental protocol
---------------------
Datasets for each subject:
- young controls: 2 drug sessions  alprazolam (Xanax 1mg) + placebo (zyrtec)
- aged subjects: 1 drug session  alprazolam (xanax 1mg)
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

------------------------------------------------------------
AVAILABLE SCRIPTS & FUNCTIONS
------------------------------------------------------------

EEG: TEPs from M1 and AG
------------------------------------------------------------
preprocessing:
- GABA_EEG_import.m
- GABA_initialize logfile.m
- GABA_EEG_parse_events.m
- GABA_EEG_preprocess1.lwscript
- GABA_EEG_preprocess2.lwscript

TEP analysis
- GABA_EEG_amplitude.m
- GABA_EEG_meanamp.m
- GABA_EEG_topoplot.m

EMG: MEPs from M1
------------------------------------------------------------
-

RS-EEG
------------------------------------------------------------
-
