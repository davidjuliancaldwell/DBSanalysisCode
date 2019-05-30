### Code to analyze evoked potentials from intraoperative DBS experiments, as well as paired pulse stimulation experiments in humans.

This repository contains MATLAB and R code to analyze the data from different sorts of DBS recordings, as well as paired pulse experiments. There exists analysis code for paired pulse experiments from both DBS and epilepsy patients.

### MATLAB extraction of evoked potentials from paired pulse experiments

The main script is ***master_script_betaStim.m***, in the ***DBS_EP_PairedPulse*** folder, which calls other sub scripts and analysis functions. This generates time series plots of the different evoked potentials, as well as CSV files for further statistical analyses.

---

### R analysis of evoked potentials

The ***R_code*** folder under ***DBS_EP_PairedPulse*** contains the R scripts required to fit linear mixed models and generate statistical plots after the data structure generated from the MATLAB code ***master_script_analyze_EP.m*** has been run. Each individual subject has various EP parameters that need to be set in ***prepare_EP_blocks.m***.

***screenBadChans*** is a boolean variable, which when enabled, allows for the visualization of each individual evoked potential to allow for manual rejection of bad trials.

---

After running the MATLAB analysis, the R analysis can be run.  

***R_data*** should contain the data files written out by the above MATLAB code.

***R_config_files*** contains the individual subject configuration files for the R analysis scripts.

***dose_response_R_script_trim_conditions.R*** produces the analysis contained within my thesis across all the subjects

***baseline_variability_3d413.R*** produces the results for the patient of baseline variability under differing levels of anesthesia.

***length_conditioning_a23ed.R*** analyzes the results from the different length conditioning sessions.

---
### Intraoperative Paired Pulse Experiments

The folder ***DBS_EP_PairedPulse/intraoperative*** contains code to generate time series evoked potential plots during the running of experiments. ***intraop_DBS_patients*** is for DBS patients, while ***intraop_epilepsy_patients*** is for epilepsy patients.

The MATLAB script ***intraop_compare_EP_master_script.m*** needs to modified for each particular subject.  

---

David J Caldwell, BSD-3 License
