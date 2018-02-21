# smFISH-analysis
Analyze smFISH data

Original scripts (of Stefan):
* completefit_differentiation_290118.m
* dotfitting_parameters.m
* fftfilter_nowrap_fast_pad.m
* findcandidates_3d_manual.m
* improcess_2015.m

<< Dependencies >>
Main file: completefit_differentiation_290118.m
1. Select files (figures, stacks)
2. Choose dot fitting parameters [dotfitting_parameters.m]
3. Process images [improcess_2015.m]
    * Background subtraction [fftfilter_nowrap_fast_pad.m]
    * Denoising [fftfilter_nowrap_fast_pad.m]
4. Find candidates (smFISH dots) and select the right threshold manually [findcandidates_3d_manual.m] 
    * Project 3D image stack to a 2D image using maximum projection
    * Calculate noise from maximum projection
    * Thresholding > find candidates
    * Measure properties of each candidate
    * Save candidates on image
5. (Dependent script [assignpeakstocells_2015] not included: assign dots to cells using DAPI image)

Script 'smFISH_analysis_esmee_v1': all necessary parts of the scripts are written here in one script,
so you can easily go through all the steps yourself.
