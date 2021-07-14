clear all
%
% t-CycIF ASHLAR ANALYSIS PIPELINE 
% Before starting: Stitch and register the fields into an "ome.tif" by using ASHLAR (https://doi.org/10.1101/2021.04.20.440625)
% Step 1:(Step1_fields_preilastik.m) Cuts the "ome.tif" into fields of a size specified by the user and cut crops out of each field 
% for Ilastik training. Create full stacks of all cycles and channels. 
% Step 2: (Step2b_filtercrops_v2.m) OPTIONAL Omits fields that have no cells. 
% Please note that if the stacks produced by Matlab are not recognized as multichannel image % the addition 
% RUN_Step2_changetifffromZtoXhannel_yescrops_O2.ijm needs to be run in Fiji/ImageJ (https://imagej.net/software/fiji/)
% Step 2: (needs to be run externally through Ilastik https://www.ilastik.org/) Create segmentation probabilities using Ilastik
% Step 3: (Step3_segmentfromilastik_v2_LabelAndDilate.m) Segment based on probabilities produced by Ilastik
% Step 4: (Step4_CycIF_measurements_ilastik_v2.m) Makes measurements of signal and foci segmentation 
% 
% In order to begin the analysis a set of parameters are needed to be
% defined by the user 
% 1) where the files are located and how the names are formatted
% 2) parameters to choose the size of the field and number of crops per field
% 3) parameters for cell segmentation - e.g. cell size, probability thresholds

%
%%% IMPORTANT: Ilastik Output Filename Format for Export Options
% {dataset_dir}/Ilastik_Probabilities/{nickname}_Probabilities.tif

% 1) THE MASTER OME.TIF DATA FOLDER
filename = [];
filename.folders.main = '\CCR-001-TCF1CCR6\';

% 2) SLIDE SPECIFIC PARAMETERS:
filename.tissues = {'CASE1','CASE7','CASE8','CASE9','CASE10','CASe11','CASe12','CASE13','CASE14','CASE15','CASE16','CASE17','CASE18','CASE19','CASE20','CASE21'};
                %Name of Ashlared image without '.ome.tif' ending 
            
filename.cycles = 9; % # of cycles to analyse
filename.maxround = filename.cycles*4;
filename.ilastikround = 36; % CyCles 1 - filename.ilastikround are the rounds that will be used for ilastik

% 3) USER DESIRED PARAMETERS 
filename.sizefield = 5000; %Size of field desired for a square 
filename.crops = 1; %# of cropped fields desired per field 
options.DAPIbkgd = 100; % thresh for DAPI used to check if we keep the field
options.DAPIbkgd_crops75prctile = options.DAPIbkgd*5;

% 4) OPTIONS  
% Step 3: SEGMENTATION OPTIONS 
options.nuc = 1;        % Channel nucleus Ilastik probability is in
options.cyt = 2;        % Channel cytoplasm Ilastik probability is in
options.backgr = 3;     % Channel background Ilastik probability is in 
options.cellsize = 9;   % Estimated radius of cells in pixels
options.bkgdprob_min = 0.15;  % Thresh for backgorund Ilastik probability
options.cytprob_min  = 0.33;  % Thresh for backgorund Ilastik probability
options.nucprob_min  = 0.7;   % Thresh for backgorund Ilastik probability
options.max_prob = 65535;     % Maximum Ilastik probability for rescaling
options.pausetime = 0;

% Step 4: MEASUREMENTS AND FOCI OPTIONS 
options.date = '20201106';  % add the date here

% Step 5 which slice to use for the montage for ROI
options.pyramidlevel = 3;
options.chan_pyramid = 12; 

% 5) FOLDER OUTPUT PARAMETERS: DO NOT EDIT 
filename.folders.ashlared = '\DATA\';   
filename.ashlarsubfol = '\registration\';
filename.ashlarsuffix = '.ome.tif';
filename.folders.output = '\ANALYSIS\'; 
filename.folders.fullstacks = '\FullStacks\';
filename.folders.coordinates = '\Coordinates_FullStacks\';

filename.folders.ilastikstacks = '\IlastikStacks\';
filename.folders.cropfol = '\CroppedData\';
filename.folders.ilastikprob= '\IlastikProb\';
filename.folders.ilastikseg = '\IlastikSeg\';
filename.ilastiksuffix = '_Probabilities.tif';

filename.folders.ometiff = '\DATA\';
filename.folders.montagelowres = '\MontageforROI_Lv3\';

filename.folders.results = '\RESULTS\';
filename.folders.ROI = '\ROIs\'; 
filename.suffix = '.tif';
% filename.folders.montage = 'MontageforROI\'; 
filename.dim = ['%02d']; %Delimeter 
filename.overwrite_seg = 0;

%% MAKE SURE TO COMMENT OUT THE OTHER STEPS THAT YOU AREN'T AT YET!!
% Runs Step 1: (Creating FullStacks, CroppedStacks, Montages) 
Step1_fields_preilastik(filename,options) 
%% clean up the crops to avoid emptry crops in Ilastik training
Step2b_filtercrops_v2(filename,options) 
% after this step need to perform Ilastik step to create probability maps
%% Runs Step 3 (Single cell segmentation)
Step3_segmentfromilastik_v2_LabelAndDilate(filename, options) 
%% Runs Step 4 (Single cell measurments)
Step4_CycIF_measurements_ilastik_v2(filename, options) 
%
