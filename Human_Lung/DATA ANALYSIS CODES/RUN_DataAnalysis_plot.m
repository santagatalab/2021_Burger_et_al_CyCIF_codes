%%
clear all
codedir = 'Y:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\CycIF Codes\Utility Functions';
addpath(codedir)
basefolder = 'Y:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\2020-10_Human_Lung_DataAnalysis\2020-07_HumanLung_Tcellpheno\';
analfolder = [basefolder 'ANALYSIS\'];
resufolder = 'Results_20201004\'; 
date = '20201004';

load([analfolder resufolder 'Results_Morp_' date '.mat'])
% load([analfolder resufolder 'Results_ROI_' date '.mat'])
load([analfolder resufolder 'Results_Norm_' date '.mat']);
load([analfolder resufolder 'Results_Filt_' date '.mat']);
load([analfolder resufolder 'Results_Settings_' date '.mat'])
load([analfolder resufolder 'Results_CellType_' date '.mat'])
load([analfolder resufolder 'Results_Nets_' date '_dist30.mat'])
options.figOpt = 0;

filename.basefolder = basefolder;
filename.analfolder = analfolder;
filename.resufolder = resufolder; 
options.date = date;

save([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'],'filename','options')

forplot.grouplist = {'1','2','3','4','5','6'};
forplot.group_mat = 1:6;

%%

type = [1111, 1112, 1113];
chans = [14 12 16];
mat = [NormResults.MedianNucNorm(:,chans(1:2)), NormResults.MedianCytNorm(:,chans(3))];

figure
for i = 1:length(type)
    for tis = 1:6
        ind_type = MorpResults.Indexes == tis & CellType.Matrix(:,4) == type(i);
        mat_rand = datasample(mat(ind_type,:),1000);
        subplot(3,6,6*(i-1)+tis)
        scatter(mat_rand(:,1),mat_rand(:,2),3,'filled','MarkerFaceAlpha',0.9) % mat(ind_type,3)
        axis([-2000 2000 -2000 2000])
        caxis([-1000 1000])
        colormap(NormResults.colorMap)
        title(CellType.names(find(CellType.codes==type(i))))
        xlabel(options.Markers{chans(1)})
        ylabel(options.Markers{chans(2)})
    end
end
%%
figure
for i = 1:length(type)
    ind_type = CellType.Matrix(:,4) == type(i);
    ksdensity(NormResults.MedianCytNorm(ind_type,find(strcmp(options.Markers,'GzmB'))))
    hold on
    
end

%% distribution of lymphocytes in tissues

celltypes = {'B','T reg','T helper','T cytotox','Epithelial','Other'};
mat = [];
for g = 1:size(forplot.group_mat,2)
    index = MorpResults.Indexes == forplot.group_mat(g) & Filter.all(:,1) == 1 ;
    totcells = sum(index);
    for type = 1:length(celltypes)
        code = CellType.codes(strcmp(celltypes{type},CellType.names));
        ind_type = CellType.Matrix(:, ceil(log10(double(code+1)))) == code;
        mat(g,type)=sum(ind_type & index)/totcells;
    end
end

bar(mat,'stacked')
legend('B','T reg','T helper','T cytotox','Epithelial','Other')

%% characterization of CD8 cells
celltypes = {'T cytotox','Epithelial'};
% total CD8 cells per tissue
figure
clear vect
for type = 1:length(celltypes)  
    for g = 1:size(forplot.group_mat,2)
        index = MorpResults.Indexes == forplot.group_mat(g) & Filter.all(:,1) == 1 ;
        totcells(g) = sum(index);
        
        
        code = CellType.codes(strcmp(celltypes{type},CellType.names));
        ind_type = CellType.Matrix(:, ceil(log10(double(code+1)))) == code;
        vect(type,g) = sum(ind_type & index);
    end
end
subplot(1,2,1)
bar(vect(1,:)./totcells)
title('T cyto norm by total cell #')
subplot(1,2,2)
bar(vect(1,:)./vect(2,:))
title('T cyto norm by epithelial cell #')

tcf = strcmp('TCF1-7',options.Markers);
ccr = strcmp('CCR6',options.Markers);
clear mat
for g = 1:size(forplot.group_mat,2)
    index = MorpResults.Indexes == forplot.group_mat(g) & ...
            prod(Filter.all(:,max(ccr,tcf)),2) == 1 & ...
            CellType.Matrix(:,4) == 1113 ;
    TCF1 = NormResults.MedianNucNorm(index,tcf);
    CCR6 = NormResults.MedianCytNorm(index,ccr);
    mat(g,:) = [sum(TCF1<0 & CCR6 < 0) sum(TCF1<0 & CCR6 >= 0) ...
              sum(TCF1>=0 & CCR6 < 0) sum(TCF1>=0 & CCR6 >= 0)  ]/length(TCF1);
    DoublePosEst(g) = (sum(TCF1>=0)/length(TCF1))*(sum(TCF1>=0)/length(CCR6));
end
figure
subplot(1,2,1)
bar(mat,'stacked')
subplot(1,2,2)
scatter(DoublePosEst,mat(:,4),25,'filled')
hold on
plot([0 0.2],[0 0.2],'--k')

%%
clear channels
mark = {'PD-1','GzmB','CD8a '};
for m = 1:length(mark)
    channels(m) = find(strcmp(options.Markers,mark{m}));
end
for i = 1:2
    index = MorpResults.Indexes == 1 & ...
              prod(Filter.all(:,channels),2) == 1 & ...
              CellType.Matrix(:,4) == 1113;
    if i == 1
        index = index & NormResults.MedianCytNorm(:,find(strcmp(options.Markers,'GzmB'))) > 500;
    elseif i == 2
        index = index & NormResults.MedianCytNorm(:,find(strcmp(options.Markers,'PD-1'))) > 1000;
    end
    
    cen_cd8 = [MorpResults.X(index) MorpResults.Y(index) ];
    sum(index)

    cen_cd8 = datasample(cen_cd8,32);

    tilesize = 5000;
    fullstackloc = 'Y:\sorger\data\RareCyte\Claire\2020_07_Human_Lung_CCR6_TCF1\ANALYSIS\FullStacks\CASE1';
    cropsize = 20;
    % channels = [7 8 12];

    RetrieveSingleCells(cen_cd8,tilesize,channels,fullstackloc,cropsize)
end


%% DENSITY OF CELLS

bin_size  = 500;
min_cell = 50;

th_nets_cellnum = 25;

% loop around each mouse
for g = 1:size(forplot.group_mat,2)
    index = MorpResults.Indexes == forplot.group_mat(g) ;
    cent_x = double(ceil(MorpResults.X(index)/bin_size));
    cent_y = double(ceil(MorpResults.Y(index)/bin_size));

    edges{1} = linspace(0,max(cent_x),max(cent_x)+1)+0.5;
    edges{2} = linspace(0,max(cent_y),max(cent_y)+1)+0.5;
    [N,c] = hist3([cent_x cent_y],'Edges',edges);
    clear D
    for type = 1:max(CellType.index)
        index_type = CellType.Matrix(:,CellType.layer(type))==CellType.codes(type) & index;
        cent_x_type = double(ceil(MorpResults.X(index_type)/bin_size));
        cent_y_type = double(ceil(MorpResults.Y(index_type)/bin_size));
        [temp,c] = hist3([cent_x_type cent_y_type],'Edges',edges);
        D{type} = temp(N(:)>min_cell);
        D_mean{type}(1,g) = mean(D{type}); 

    end
       
end

figure(11)
for type = 1:max(CellType.index)
    subplot(4,4,type)
    mouse_stat_plot(D_mean{type})
    set(gca,'XTick',1:4,'XTickLabel',forplot.grouplist,'XTickLabelRotation',45)
    title(CellType.names{type})
    ylabel('cells density')
end
set(gcf,'Units','normalized','Position',[0.1 0.1 0.3 0.8 ])
pause(3)
suptitle('Whole Lung Area')