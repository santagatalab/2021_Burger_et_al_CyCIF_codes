clear all
%%%
filename.basefolder = '2020-10_HumanLung\';
filename.suffix = '_Results_20201114.mat';
filename.analfolder = [filename.basefolder 'ANALYSIS\'];
filename.resufolder = 'Results_20201114\'; 

filename.folders = {'CASE1', 'CASE7','CASE8','CASE9','CASE10', ...
                    'CASE11','CASE12','CASE13','CASE14','CASE15',...
                    'CASE16','CASE17','CASE18','CASE19','CASE20','CASE21'}; 
                
filename.tissues = filename.folders;

options.Markers =  {  'DAPI0','CD4','CCR6','Granzyme B', ...
                      'DAPI1','TCF1','FOXP3','CD8a', ...
                      'DAPI2','TTF1','PD-L1','CD20',...
                      'DAPI3','TIM-3','CD45','PD-1',...
                      'DAPI4','CD163','CD68','Ki-67',...
                      'DAPI5','MHC-II','CD3D','MHC-I', ...
                      'DAPI6','PCNA','ASMA','Vimentin',...
                      'DAPI7','CD16','Keratin','CD-14',...
                      'DAPI8','CD19','CCR6','CD103'}; 
                 
options.maxround = length(options.Markers);
options.magnification = 20;
options.FigOpt = 0;
options.date = '20201202';

% STEP 2: additional parameters for filtering          
options.Filtering.folder = [filename.analfolder filename.resufolder 'Step2_'];
options.Filtering.Index_Names = filename.tissues;
options.Filtering.thresholds.foldDAPI_th = 0.5;
options.Filtering.thresholds.absDAPI_th = 9;
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
      
save([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'],'filename','options')

%% Step 1 aggregate tissues together
clearvars -except filename options
[AggrResults, MorpResults] = PreProcess_Step1_Aggregation_v2(filename, options, options.FigOpt);

save([filename.analfolder filename.resufolder 'Results_Aggr_' options.date '.mat'],'AggrResults','-v7.3');
save([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'],'MorpResults');

%% Step 2 - Filter the data

clearvars -except filename options
load([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Aggr_' options.date '.mat'])
options.FigOpt = 0;

[Filter,report] = PreProcess_Step2_Filter(AggrResults, options.Filtering, options.FigOpt);
save([filename.analfolder filename.resufolder 'Results_Filt_' options.date '.mat'],'Filter','report')

%%
% close all
clearvars -except filename options
load([filename.analfolder filename.resufolder 'Results_Filt_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Aggr_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'])

rng(11)
close all

addpath([filename.analfolder 'norm']);
    % ^ where the normalization codes are found
options.FigOpt = 1;


%%% choose channels
bad_channels = [3]; % channels not to normalize (result column will be 0s)
nondapi_channels = setdiff(1:length(options.Markers), 1:4:length(options.Markers));
channels = setdiff(nondapi_channels,bad_channels);
cyt_channels = [];

% set up normalization
priors = zeros(length(options.Markers),1)+1;

    overexpr = [0    0    0    0        ...  'DAPI0','CD4','CCR6','Granzyme B', ... 
                0    0.25 0    0        ...  'DAPI1','TCF1','FOXP3','CD8a', ...
                0    0    0    0        ...  'DAPI2','TTF1','PD-L1','CD20',...
                0    0    0    0        ...  'DAPI3','TIM-3','CD45','PD-1',...
                0    0.17 0.1  0        ...  'DAPI4','CD163','CD68','Ki-67',...
                0    0    0    0        ...  'DAPI5','MHC-II','CD3D','MHC-I', ...
                0    0    0    0        ...  'DAPI6','PCNA','ASMA','Vimentin',...
                0    0    0    0        ...  'DAPI7','CD16','Keratin','CD-14',...
                0    0    0    0        ...  'DAPI8','CD19','CCR6','CD103'
                ]; 
          
            
options.Norm.Priors = priors;
options.Norm.OverExpr = overexpr;
options.Norm.CellNum = min([50000, length(AggrResults.MedianNucSign)]);
options.Norm.Channels = channels;

filter = Filter.all;
data_nuc = log2(double(AggrResults.MedianNucSign)+1);
data_cyt = log2(double(AggrResults.MedianCytSign)+1);

% normalize Cytoplasmic signal
options.Norm.IsNuc = 0;
[NormResults.MedianCytNorm, cutoffs_cyt, mults_cyt] = norm_main(data_cyt, filter, options.Markers, options.Norm);
options.Norm.OverExpr_Cyt = options.Norm.OverExpr;

% normalize Nuclear signal

% adjust CD4 nuclear normalization
options.Norm.OverExpr(2) = 0.25;
options.Norm.OverExpr(12) = 0.1;
options.Norm.OverExpr_Nuc = options.Norm.OverExpr;

options.Norm.IsNuc = 1;
[NormResults.MedianNucNorm, cutoffs_nuc, mults_nuc] = norm_main(data_nuc, filter, options.Markers, options.Norm);

close all

%% Normalized CCR6 by tissue
CCR6r = 3;
options_adhoc = options.Norm;
options_adhoc.Channels = CCR6r;
options_adhoc.OverExpr(3) = 0.15;
filter = Filter.all;
data_for_norm = zeros(size(AggrResults.MedianNucSign))+NaN;
figure
for t = 1:16
    % select cells from single tissue
    index = MorpResults.Indexes == t;
    filt = Filter.all(:,CCR6r) == 1;
    filter_temp = Filter.all(index & filt,:);
    data_temp = log2(double(AggrResults.MedianNucSign(index & filt,:))+1);
    [n,h] = hist(data_temp(:,CCR6r),linspace(1,16,200));
    [~,p] = max(n);
    data_for_norm(index,CCR6r) = log2(double(AggrResults.MedianNucSign(index,CCR6r)+1)) - h(p);
    
    [n,h] = hist(data_temp(:,CCR6r)-h(p),linspace(1,16,200)-h(p));
    plot(h,n/sum(n))
    hold on
    title('CCR6 nuclear signals normalized by negative distribution')
end

options_adhoc.IsNuc = 1;
[data_norm, ~, ~] = norm_main(data_for_norm, filter, options.Markers, options_adhoc);
NormResults.MedianNucNorm(:,CCR6r) = data_norm(:,CCR6r);

options_adhoc.OverExpr(3) = 0.14;

data_for_norm = zeros(size(AggrResults.MedianNucSign))+NaN;
figure
for t = 1:16
    % select cells from single tissue
    index = MorpResults.Indexes == t;
    filt = Filter.all(:,CCR6r) == 1;
    filter_temp = Filter.all(index & filt,:);
    data_temp = log2(double(AggrResults.MedianNucSign(index & filt,:))+1);
    [n,h] = hist(data_temp(:,CCR6r),linspace(1,16,200));
    [~,p] = max(n);
    data_for_norm(index,CCR6r) = log2(double(AggrResults.MedianCytSign(index,CCR6r)+1)) - h(p);
    
    [n,h] = hist(data_temp(:,CCR6r)-h(p),linspace(1,16,200)-h(p));
    plot(h,n/sum(n))
    hold on
    title('CCR6 cytoplasmic signal per tissue normalized by negative distribution')
end

options_adhoc.IsNuc = 0;
[data_norm, ~, ~] = norm_main(data_for_norm, filter, options.Markers, options_adhoc);
NormResults.MedianCytNorm(:,CCR6r) = data_norm(:,CCR6r);

%% save
redmap = [linspace(0,255,128) zeros(1,128)+255 ]/255;
blumap = [zeros(1,128)+255 flip(linspace(0,255,128))]/255;
gremap = [linspace(128,255,128) flip(linspace(128,255,128))]/255;
NormResults.colorMap = [redmap' gremap' blumap'];
NormResults.CellID = (1:size(NormResults.MedianNucNorm,1))';

NormResults.MedianNucNorm = int16(round(NormResults.MedianNucNorm*1000,0));
NormResults.MedianCytNorm = int16(round(NormResults.MedianCytNorm*1000,0));
NormResults.nuc_add_fact = int16(1000.*cutoffs_nuc); 
NormResults.nuc_mult_fact = int16(1000.*mults_nuc);
NormResults.cyt_add_fact = int16(1000.*cutoffs_cyt); 
NormResults.cyt_mult_fact = int16(1000.*mults_cyt);

save([filename.analfolder filename.resufolder 'Results_Norm_' options.date '.mat'],'NormResults','options','filename');
save([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'],'filename','options');

