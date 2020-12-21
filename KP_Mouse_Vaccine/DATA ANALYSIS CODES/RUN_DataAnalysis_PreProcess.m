clear all
%%%
codedir = 'Y:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\CycIF Codes\Utility Functions';
addpath(codedir)
% addpath('...\normalization\functions\')
    % ^ where the normalization codes are found
%%%
filename.basefolder = 'Y:\sorger\data\RareCyte\Claire\CCR-002-ImmuneVaccine\';
filename.suffix = '_Results_20201117.mat';
filename.analfolder = [filename.basefolder 'ANALYSIS\'];
filename.resufolder = 'Analysis_Results_20201105\'; 
filename.roifolder  = 'ROIs\';
filename.montfolder = 'MontageforROI_Lv3\';
filename.montsuffix = '_montage.tif';

filename.folders =  {'MB_1117','MB_1118','MB_1119','MB_1120','MB_1121' ...
                    ,'MB_1122','MB_1123','MB_1124','MB_1125','MB_1126' ...
                    ,'MB_1127','MB_1128','MB_1129','MB_1130','MB_1131' ...
                    };
                
filename.tissues = filename.folders;              

for i = 1:length(filename.folders)
    options.MouseNum(i) = str2num(filename.folders{i}(4:7));
    options.MouseGroup(i) = ceil(double(options.MouseNum(i))/1123);
end

options.Markers =  {'DAPI 0',   'Bkgd488', 'Bkgd555',   'Bkgd647', ...
                    'DAPI 1',   'TTF1',    'CD45R-B220','CD45', ...
                    'DAPI 2',   'FOXP3',   'CD4',       'CD8a',...
                    'DAPI 3',   'GzmB',    'Casp3' ,    'Ki67'}; 
                  
options.maxround = length(options.Markers);
options.magnification = 20;
options.FigOpt = 0;
options.date = '20201119';

% STEP 2: additional parameters for filtering          
options.Filtering.folder = [filename.analfolder filename.resufolder 'Step2_'];
options.Filtering.Index_Names = filename.tissues;
options.Filtering.thresholds.foldDAPI_th = 0.5;
options.Filtering.thresholds.absDAPI_th = 12;
options.Filtering.thresholds.solidity = 0.8;
options.Filtering.thresholds.area_low = 50;    % was 50   in previous analysis
options.Filtering.thresholds.area_high = 1000;  % was 2000 in previous analysis..
options.Filtering.maxround = options.maxround;

% STEP 3: additional parameters for normalization
options.Norm.Reps = 3;
options.Norm.FigSettings.FigFlag = 1; % to save
options.Norm.FigSettings.Folder = [filename.analfolder filename.resufolder 'Step3_NormPrints\' ];
options.Norm.FigSettings.Debug = 1; % to view
options.Norm.Channels = setdiff(1:length(options.Markers), 1:4:length(options.Markers));  % non-dapi channels
% default to zeros
options.Norm.Priors = zeros(length(options.Markers),1); 
options.Norm.OverExpr = zeros(length(options.Markers),1);
options.Norm.CellNum = 50000;

% % STEP 4: roi extracting
% options.tifnames.Aiforia = {'TumorGrade1.tif','TumorGrade2.tif','TumorGrade3.tif'};
options.PathLabels = {'Tumors','Spleen'};
options.ROI.pyramidlevel = 3;
% 
% % STEP 5: ROI Masking
% options.ROI.pyramidlevel = 2; % where full size is 0
% options.ROI.FigOpt = true;
% options.ROI.date = '20200424';                 

save([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'],'filename','options')
disp('DONE')
%% Step 1 aggregate tissues together
clearvars -except filename options
[AggrResults, MorpResults] = PreProcess_Step1_Aggregation_v2(filename, options, options.FigOpt);

save([filename.analfolder filename.resufolder 'Results_Aggr_' options.date '.mat'],'AggrResults','-v7.3');
save([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'],'MorpResults');
disp('Aggregation DONE')
%% Step 2 - Filter the data

clearvars -except filename options
load([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Aggr_' options.date '.mat'])
options.FigOpt = 0;

[Filter,report] = PreProcess_Step2_Filter(AggrResults, options.Filtering, options.FigOpt);
save([filename.analfolder filename.resufolder 'Results_Filt_' options.date '.mat'],'Filter','report')
disp('Filtering DONE')
%%
% close all
clearvars -except filename options
load([filename.analfolder filename.resufolder 'Results_Filt_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Aggr_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'])

rng(11)
close all
addpath([filename.analfolder 'DATA ANALYSIS CODES\norm']);
    % ^ where the normalization codes are found
options.FigOpt = 1;


%%% choose channels
bad_channels = []; % channels not to normalize (result column will be 0s)
nondapi_channels = setdiff(1:length(options.Markers), 1:4:length(options.Markers));
channels = setdiff(nondapi_channels,bad_channels);
cyt_channels = [];

close all

% set up normalization
priors = zeros(length(options.Markers),1); 
    priors([ 6 7 8 11 12 14 15 16 ]) = +1;
    overexpr = [0   0  0   0        ...  'DAPI 0',   'TTF1',    'B2m' ,    'pStat1', ...
                0   0  0   0.28        ...  'DAPI 1',   'FOXP3',    'CD4' ,    'CD8a', ...
                0   0.25  0   0        ...  'DAPI 2',   'GzmB',     'Perforin','Ki67',...
                0   0.12  0   0        ...
                ]; 
     
                   
options.Norm.Priors = priors;
options.Norm.OverExpr = overexpr;
options.Norm.CellNum = min([50000, length(AggrResults.MedianNucSign)]);
options.Norm.Channels = channels;
filter = Filter.all;
data_nuc = log2(double(AggrResults.MedianNucSign)+1);
data_cyt = log2(double(AggrResults.MedianCytSign)+1);

% normalize
options.Norm.IsNuc = 1;
[NormResults.MedianNucNorm, cutoffs_nuc, mults_nuc] = norm_main(data_nuc, filter, options.Markers, options.Norm, options.Norm.FigSettings);

options.Norm.IsNuc = 0;
[NormResults.MedianCytNorm, cutoffs_cyt, mults_cyt] = norm_main(data_cyt, filter, options.Markers, options.Norm, options.Norm.FigSettings);

%%% save
redmap = [linspace(0,255,128) zeros(1,128)+255 ]/255;
blumap = [zeros(1,128)+255 flip(linspace(0,255,128))]/255;
gremap = [linspace(128,255,128) flip(linspace(128,255,128))]/255;
NormResults.colorMap = [redmap' gremap' blumap'];
NormResults.CellID = uint16((1:size(NormResults.MedianNucNorm,1))');

NormResults.MedianNucNorm = int16(round(NormResults.MedianNucNorm*1000,0));
NormResults.MedianCytNorm = int16(round(NormResults.MedianCytNorm*1000,0));
NormResults.nuc_add_fact = int16(1000.*cutoffs_nuc); 
NormResults.nuc_mult_fact = int16(1000.*mults_nuc);
NormResults.cyt_add_fact = int16(1000.*cutoffs_cyt); 
NormResults.cyt_mult_fact = int16(1000.*mults_cyt);

save([filename.analfolder filename.resufolder 'Results_Norm_' options.date '.mat'],'NormResults');
save([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'],'filename','options');

close all
for c = 1:length(channels)
    figure
    subplot(1,2,1)
    [n,h]=ksdensity(NormResults.MedianNucNorm(Filter.all(:,channels(c))==1,channels(c)));
    plot(h,n)
    hold on
    plot([0 0],[0 max(n)])
    xlim([-2000 2000])
    subplot(1,2,2)
    [n,h]=ksdensity(NormResults.MedianCytNorm(Filter.all(:,channels(c))==1,channels(c)));
    plot(h,n)
    hold on
    plot([0 0],[0 max(n)])
    xlim([-2000 2000])
end

%% Step 4: use ROI montages to extract annotations from other sources
% this step depends on where the masks come from and will need to output a
% structure with all the mask types (ie tumor, blood, spleen) containing
% binary or bwlabeled masks that fit perfectly over the original tissue
clearvars -except filename options

PreProcess_Step4_ExtractROIs(filename,options)

%% Step 5: find whether cells belong to ROIs
clearvars -except filename options
load([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'])

ROIResults.TumorIndex = zeros(size(MorpResults.X));
ROIResults.SpleenIndex = zeros(size(MorpResults.X));

for i = 1:length(filename.folders)
    i
    % load the images
    tumor_file = [filename.analfolder filename.roifolder filename.folders{i} '_' options.PathLabels{1} '.tif'];
    spleen_file = [filename.analfolder filename.roifolder filename.folders{i} '_' options.PathLabels{2} '.tif'];
    ROI.Tumors{i} = uint16(imread(tumor_file));
    ROI.Spleen{i} = uint16(imread(spleen_file));
    
    % first find all the cells that belong to the ROI
    index = MorpResults.Indexes == i;
    % get the centroids of the cells, rescale and reshape them
    X_tissue = MorpResults.X(index);
    Y_tissue = MorpResults.Y(index);
    % rescale
    X_tissue = round(double(X_tissue)/(2^(options.ROI.pyramidlevel-1)),0);
    Y_tissue = round(double(Y_tissue)/(2^(options.ROI.pyramidlevel-1)),0);
    % make sure there are no zeros
    X_tissue(X_tissue==0) = 1;
    Y_tissue(Y_tissue==0) = 1;
    % find indices
    xy_index = sub2ind(size(ROI.Tumors{i}),Y_tissue,X_tissue);    
       
    ROIResults.TumorIndex(index)  = ROI.Tumors{i}(xy_index) ;    
    ROIResults.SpleenIndex(index) = ROI.Spleen{i}(xy_index) ;

    figure
    subplot(1,2,1)
    scatter(X_tissue,Y_tissue,5,ROIResults.SpleenIndex(index),'filled')
    subplot(1,2,2)
    scatter(X_tissue,Y_tissue,5,ROIResults.TumorIndex(index),'filled')
end

save([filename.analfolder filename.resufolder 'Results_ROI_' options.date '.mat'],'ROIResults')
save([filename.analfolder filename.resufolder 'Images_ROI_' options.date '.mat'],'ROI','-v7.3')

%%
% % 
% % for i = 1:length(filename.folders)
% %     ROI.Aiforia{i} = [];
% %     for grade = 1:length(options.tifnames.Aiforia)
% %         Grade{grade} = uint16(imread([filename.analfolder filename.roifolder filename.folders{i}(5:end) '_' options.tifnames.Aiforia{grade}]));
% %       
% %         if sum(ROI.Aiforia{i}) == 0
% %             ROI.Aiforia{i} = Grade{grade}*grade;
% %         else
% %             ROI.Aiforia{i} = ROI.Aiforia{i} + Grade{grade}*grade;
% %         end
% %     end
% %     if options.FigOpt == 1
% %         figure
% %         imshow(ROI.Aiforia{i},[])
% %         title(filename.folders{i}(5:end))
% %     end
% %     ROI.Tumors{i} = [];
% %     ROI.Spleen{i} = [];
% %     ROI.Blood{i} = [];
% %     for sc = 1:length(options.PathLabels)
% %         Score{sc} = int16(imread([filename.analfolder filename.roifolder filename.folders{i}(5:end) '_' options.PathLabels{sc}]));
% %         if strcmp(options.PathLabels{sc}(1:end-4),'Tumors')
% %             ROI.Tumors{i} = uint16(Score{sc});
% %         elseif strcmp(options.PathLabels{sc}(1:end-4),'Blood')
% %             ROI.Blood{i} = uint16(Score{sc});
% %         elseif strcmp(options.PathLabels{sc}(1:end-4),'Spleen')
% %             ROI.Spleen{i} = uint16(Score{sc});
% %         end
% %     end
% %     if options.FigOpt == 1
% %         figure
% %         imshow(ROI.Tumors{i}-ROI.Blood{i}-ROI.Spleen{i}*2,[-2 10])
% %         title(filename.folders{i}(5:end))
% %     end
% % end
% % save([filename.analfolder filename.resufolder 'Results_ROI_' options.date '.mat'],'ROI')
% % % 
% % %% Step 5: find whether cells belong to ROIs
% % clearvars -except filename options
% % load([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'])
% % load([filename.analfolder filename.resufolder 'Results_ROI_' options.date '.mat'])
% % load([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'])
% % 
% % MorpInput.X = MorpResults.X; MorpInput.Y = MorpResults.Y;
% % MorpInput.Indexes = MorpResults.Indexes; clearvars MorpResults
% % 
% % ROIResults = PreProcess_Step5_ROIMasking(filename, options.ROI, ROI, MorpInput);
% % 
% % save([filename.analfolder filename.resufolder 'Results_ROI_' options.date '.mat'],'ROIResults','-append')













