%%
clear all
codedir = 'Y:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\CycIF Codes\Utility Functions';
addpath(codedir)
basefolder = 'Y:\sorger\data\RareCyte\Claire\CCR-002-ImmuneVaccine\';
analfolder = [basefolder 'ANALYSIS\'];
resufolder = 'Analysis_Results_20201105\'; 
date = '20201119';

load([analfolder resufolder 'Results_Morp_' date '.mat'])
load([analfolder resufolder 'Results_ROI_' date '.mat'])
load([analfolder resufolder 'Results_Norm_' date '.mat']);
load([analfolder resufolder 'Results_Filt_' date '.mat']);
load([analfolder resufolder 'Results_Settings_' date '.mat'])
load([analfolder resufolder 'Results_CellType_' date '.mat'])
options.figOpt = 0;

filename.basefolder = basefolder;
filename.analfolder = analfolder;
filename.resufolder = resufolder; 
options.date = date;

save([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'],'filename','options')

forplot.grouplist = {'ctrl','vax'};
forplot.group_mat = [1117:1123 0;  1124:1131];%

c = double(size(forplot.group_mat,2));
r = double(size(forplot.group_mat,1));
forplot.plotvect = repmat((1:r),c,1) + repmat( ((1:c)'-(c+1)/2)*0.05,1,r);
         
%% Compare Lymphocytes between treatments at tumor level
clear cellprop_all
celltypecode = [111 112 1111 1112 1113];
cellprop_all = [];
for i = 1:length(celltypecode)
    celltypename{i} = CellType.names{CellType.codes==celltypecode(i)};
end

% select two conditions to test
test_group1 = 1;
test_group2 = 2;

% create vector to plot cd8 frequency
cd8infilt = zeros(size(CellType.Matrix,1),1)-0.1;

for t = 1:length(celltypecode)
    cd_mat = [];
    typecode = celltypecode(t); 
    typelayer = ceil(log10(typecode+1));
    typename = CellType.names{find(CellType.codes==typecode)};

    % loop around the mouse cohort/conditions
    cellprop = zeros(size(forplot.group_mat));
    count = 0;
   
    for g = 1:size(forplot.group_mat,1)
        % loop around the mouse within each condition
        for i = 1:size(forplot.group_mat,2)
            % isolate cells in tumor areas of specific mouse
            rois_ind = find(options.MouseNum==forplot.group_mat(g,i));
            index = ismember(MorpResults.Indexes,rois_ind) & ROIResults.SpleenIndex == 0 & ROIResults.TumorIndex > 0;
            index_notum = ismember(MorpResults.Indexes,rois_ind) & ROIResults.SpleenIndex == 0;
            
            
            if sum(index) > 1000
                % count cells inside the tumor areas
                cellcount = sum(CellType.Matrix(index,typelayer)==typecode);
                cellprop(g,i) = cellcount/sum(index)*100;
                
                % count cells in the whole lung
                cellcount = sum(CellType.Matrix(index_notum,typelayer)==typecode);
                cellprop_all(g,i) = cellcount/sum(index_notum)*100;
                
                for tum = 1:max(ROIResults.TumorIndex(index))
                    count = count + 1;
                    tum_index = ismember(MorpResults.Indexes,rois_ind) & ROIResults.SpleenIndex == 0 & ROIResults.TumorIndex == tum;

                    cellcount = sum(CellType.Matrix(tum_index,typelayer)==typecode);
                    cd_mat(count,1:3) = [g cellcount/sum(tum_index) forplot.group_mat(g,i)];
                    
                    cd8infilt(tum_index) = cellcount/sum(tum_index);
                    
                end
            else
                cellprop(g,i) = NaN;
                cellprop_all(g,i) = NaN; 
            end
            
        end 
    end
    
    % print table
    nameforfile = ['20201130_KP_Vaccine_Tumor_Infiltr_' celltypename{t} '.csv'];
    csvwrite(nameforfile,cd_mat(:,[1 3 2]))
            
    nameforfile = ['20201130_KP_Vaccine_Mouse_' celltypename{t} 'cellsOverAllCells.csv'];
    forprint = sortrows([forplot.group_mat(:) cellprop_all(:)]);
    csvwrite(nameforfile,forprint)
    

    temp_mat = cellprop_all';
    median_forbar = nanmedian(temp_mat);
    figure(1)
    subplot(length(celltypecode),3,1+(t-1)*3)
    bar([1:size(forplot.group_mat,1)],median_forbar,0.2,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
    hold on
    scatter(forplot.plotvect(:),temp_mat(:),40,round(forplot.plotvect(:)),'filled')
    caxis([0.5 size(forplot.group_mat,1)])
    xlim([0.5 size(forplot.group_mat,1)+.5])
    ylim([0 max(temp_mat(:))*1.1])
    colormap('jet')
    title([typename ' cell per mouse'])
    ylabel('% of all cells')
    
    for i = 1:2
        labels_props{i} = [forplot.grouplist{i} ' ' num2str(round(median_forbar(i),3))];
    end

    set(gca,'XTick',1:size(forplot.group_mat,1),'XTickLabel',labels_props)%,'XTickLabelRotation',45)
    
    
    clear median_forbar
    for i = 1:2
        median_forbar(i) = nanmedian(cd_mat(cd_mat(:,1)==i,2));
    end

    plotvect = cd_mat(:,1) + (rand(size(cd_mat(:,1)))*0.4) - 0.2;

    figure(1)
    subplot(length(celltypecode),3,2+(t-1)*3)
    bar([1:size(forplot.group_mat,1)]+[0.375 -0.375],median_forbar,1,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
    hold on
    scatter(plotvect,cd_mat(:,2),4,cd_mat(:,1),'filled')
    caxis([0.5 size(forplot.group_mat,1)])
    xlim([0.5 size(forplot.group_mat,1)+.5])
    ylim([0 max(cd_mat(:,2))*1.1])
    colormap('jet')
    title([typename ' cell per tumor'])
    
    for i = 1:2
        labels_props{i} = [forplot.grouplist{i} ' ' num2str(round(median_forbar(i),3))];
    end
    
    set(gca,'XTick',1:size(forplot.group_mat,1),'XTickLabel',labels_props)%,'XTickLabelRotation',45)

    
    
    [ktest,p] = kstest2(cd_mat(cd_mat(:,1)==test_group1,2),cd_mat(cd_mat(:,1)==test_group2,2));
    
    lim_max = zeros(1,size(forplot.group_mat,1));
    figure(1)
    subplot(length(celltypecode),3,3+(t-1)*3)
    for i = 1:size(forplot.group_mat,1)
        [n,h] = ksdensity(cd_mat(cd_mat(:,1)==i,2));
        norm_h = n/sum(n)/(h(2)-h(1));
        plot(h,norm_h)
        hold on
        lim_max = [max(lim_max(1),max(h)) max(lim_max(2),max(norm_h))];
    end
    title(['KStest ' num2str(test_group1) ' vs ' num2str(test_group2) ' p=' num2str(round(p,4))])
    axis([0 lim_max(1) 0 lim_max(2)*1.05])
    
    if t == length(celltypecode)
        legend(forplot.grouplist,'location','best')
        legend('boxoff')
    end

end

figure(1)
set(gcf,'Units','normalized','position',[0.35 0.1 0.3 0.8])
%% plot the cd8 frequency to look at tissues

for g = 1:size(forplot.group_mat,1)
    % loop around the mouse within each condition
    for i = 1:size(forplot.group_mat,2)
        % isolate cells in tumor areas of specific mouse
        rois_ind = find(options.MouseNum==forplot.group_mat(g,i));
        index = ismember(MorpResults.Indexes,rois_ind); % & ROIResults.SpleenIndex == 0
        xdata = MorpResults.X(index);
        ydata = MorpResults.Y(index);
        cd8data = cd8infilt(index);
        
        figure
        scatter(xdata,ydata,5,cd8data,'filled')
        caxis([-0.05 0.2])
        colormap(flipud(parula))
    end
end

%% DENSITY OF CELLS

bin_size  = 500;
min_cell = 50;

th_nets_cellnum = 15;

test_group1 = 1;
test_group2 = 2;

% loop around each mouse
for g = 1:size(forplot.group_mat,1)
    % loop around the mouse within each condition
    for i = 1:size(forplot.group_mat,2)
        % isolate cells of specific mouse
        rois_ind = find(options.MouseNum==forplot.group_mat(g,i));
        clear D
        for r = 1:length(rois_ind)
            index = MorpResults.Indexes == rois_ind(r) & ROIResults.SpleenIndex == 0;
            
            cent_x = double(ceil(MorpResults.X(index)/bin_size));
            cent_y = double(ceil(MorpResults.Y(index)/bin_size));

            edges{1} = linspace(0,max(cent_x),max(cent_x)+1)+0.5;
            edges{2} = linspace(0,max(cent_y),max(cent_y)+1)+0.5;
            [N,c] = hist3([cent_x cent_y],'Edges',edges);

            for type = 1:max(CellType.index)
                index_type = CellType.Matrix(:,CellType.layer(type))==CellType.codes(type) & index;
                cent_x_type = double(ceil(MorpResults.X(index_type)/bin_size));
                cent_y_type = double(ceil(MorpResults.Y(index_type)/bin_size));
                [temp,c] = hist3([cent_x_type cent_y_type],'Edges',edges);

                if r == 1
%                     D{type} = temp(N(:)>min_cell)./N(N(:)>min_cell);
                    D{type} = temp(N(:)>min_cell) ;%./sum(N(:)>min_cell);
                elseif r == 2
%                     D{type} = [D{type}; temp(N(:)>min_cell)./N(N(:)>min_cell)] ;
                    D{type} = [D{type}; temp(N(:)>min_cell) ];%./sum(N(:)>min_cell)] ;
                end


                D_mean{type}(g,i) = mean(D{type}); 

            end
            if r == 1
                nets_num(g,i) = sum(LymphoNets.Summary{rois_ind(r)}.NumCells > th_nets_cellnum);
            elseif r == 2
                nets_num(g,i) = nets_num(g,i) + sum(LymphoNets.Summary{rois_ind(r)}.NumCells > th_nets_cellnum);
            end
        end
                
        if length(rois_ind) == 0
            for type = 1:max(CellType.index)
                D_mean{type}(g,i) = NaN;
            end
            nets_num(g,i) = NaN;
        end
        
    end
end

figure(11)
for type = 1:max(CellType.index)
    subplot(4,4,type)
    mouse_stat_plot(D_mean{type})
    set(gca,'XTick',1:size(forplot.group_mat,1),'XTickLabel',forplot.grouplist)%,'XTickLabelRotation',45)
    [~,p]=ttest2(D_mean{type}(1,:),D_mean{type}(2,:));
    title([CellType.names{type} ' p=' num2str(round(p,5))])
    ylabel('cells density')
end
set(gcf,'Units','normalized','Position',[0.1 0.1 0.3 0.8 ])
pause(3)
suptitle('Whole Lung Area')

         
%% GzmB exploration
probe_round = find(strcmp(options.Markers,'GzmB'));

% find cd8 cells in the whole mouse and in the tumors
CD8cells = CellType.Matrix(:,4) == 1113 & Filter.all(:,probe_round)>0 & ROIResults.SpleenIndex == 0;
CD8cells_intumor = CellType.Matrix(:,4) == 1113 & ROIResults.TumorIndex > 0 & Filter.all(:,probe_round)>0;

linecols = {'b','r'};
clear temp

CD8pos_total_mouse = zeros(size(forplot.group_mat)) + NaN;
CD8pos_tumor_mouse = zeros(size(forplot.group_mat)) + NaN;
CD8pos_tumor = [];

count = 0;

for g = 1:size(forplot.group_mat,1)
    gzmB_cum = [];
    tum_ind_cum = [];
    
    % loop around the mouse within each condition
    for i = 1:size(forplot.group_mat,2)
        % isolate cells of specific mouse
        rois_ind = find(options.MouseNum==forplot.group_mat(g,i));
        clear D
        temp_CD8cells = 0;
        if length(rois_ind) == 0
            gzmb_pos(g,i) = NaN;
            continue
        end
        
        cd8_mouse = sum(CD8cells & MorpResults.Indexes == rois_ind);
        cd8pos_mouse = sum(CD8cells & MorpResults.Indexes == rois_ind & NormResults.MedianNucNorm(:,probe_round)>0);
        
        cd8_tumor = sum(CD8cells_intumor & MorpResults.Indexes == rois_ind);
        cd8pos_tumor = sum(CD8cells_intumor & MorpResults.Indexes == rois_ind & NormResults.MedianNucNorm(:,probe_round)>0);
        
        CD8pos_total_mouse(g,i) = cd8pos_mouse/cd8_mouse ; 
        CD8pos_tumor_mouse(g,i) = cd8pos_tumor/cd8_tumor ;  
        
        
        % get stats per tumor
        tumor_indices = unique(ROIResults.TumorIndex(MorpResults.Indexes == rois_ind));
        
        for t = 2:length(tumor_indices)
            % find tumor cells - skip the first index as it is == 0 ie not
            % tumor lung tissue
            tum_ind = ROIResults.TumorIndex == tumor_indices(t) & MorpResults.Indexes == rois_ind ;
            % calculate num of CD8 and CD* positive for the proeb
            CD8_pertum = sum(CD8cells & tum_ind);
            CD8pos_pertum = sum(CD8cells & tum_ind & NormResults.MedianNucNorm(:,probe_round)>0);
            
            % check that there are CD8 in this tumor
            if CD8_pertum == 0
                continue
            end
            
            count = count + 1;
            CD8pos_tumor(count,1:4) = [g forplot.group_mat(g,i) tumor_indices(t) CD8pos_pertum/CD8_pertum] ;
        end
        
        gzmb_CD8cell = NormResults.MedianNucNorm(CD8cells & MorpResults.Indexes == rois_ind,probe_round);
        tumind_temp = ROIResults.TumorIndex(CD8cells & MorpResults.Indexes == rois_ind) > 0;
        
        gzmB_cum = [gzmB_cum; gzmb_CD8cell];
        tum_ind_cum = [tum_ind_cum ; tumind_temp];
        
    end
    
    figure(33)
    subplot(2,2,1)
    [n,h]=ksdensity(gzmB_cum);
    plot(h,n,linecols{g},'LineWidth',2)
    hold on
    xlim([-2500 2500])
    
    GrzBpos_prop = sum(gzmB_cum > 0)/length(gzmB_cum);
    
    if g == 1
        text(1000,(0.7)*10^(-3),'Overall GzmB+ %' ,'HorizontalAlignment','center')
    end
    text(1000,(0.7-g*0.05)*10^(-3),[forplot.grouplist{g} ' =  ' num2str(round(GrzBpos_prop*100,2)) '%'],'HorizontalAlignment','center')
    
    title(['Single cell GzmB distribution by group'])
    xlabel('Normalized GzmB CyCIF signal')
    ylabel('ks pdf (frequency)')
    
    GrzBpos_prop_loc(g,1) = sum(gzmB_cum > 0)/length(gzmB_cum);
    GrzBpos_prop_loc(g,2) = sum(gzmB_cum(tum_ind_cum > 0) > 0)/length(gzmB_cum(tum_ind_cum > 0));
    
    
end

figure(33)
subplot(2,2,1)
legend(forplot.grouplist)
% legend({[forplot.grouplist{1} ' OUT of tumors'],[forplot.grouplist{1} ' IN tumors'],[forplot.grouplist{2} ' OUT of tumors'],[forplot.grouplist{2} ' IN tumors']})

subplot(2,2,2)
bar(GrzBpos_prop_loc')
legend(forplot.grouplist)
set(gca,'XTick',1:size(forplot.group_mat,1),'XTickLabel',{'ALL','IN'})%,'XTickLabelRotation',45)

subplot(2,2,3)
mouse_stat_plot([CD8pos_total_mouse; CD8pos_tumor_mouse])
labels = {['All ' forplot.grouplist{1}],['All ' forplot.grouplist{2}],['IN ' forplot.grouplist{1}],['IN ' forplot.grouplist{2}]};
set(gca,'XTick',1:size(forplot.group_mat,1)*2,'XTickLabel',labels)
[~,p_All]=ttest2(CD8pos_total_mouse(1,:),CD8pos_total_mouse(2,:));
[~,p_IN]=ttest2(CD8pos_tumor_mouse(1,:),CD8pos_total_mouse(2,:));
title(['p-All = ' num2str(round(p_All,5)) ' & p-IN = ' num2str(round(p_IN,5))])

subplot(2,2,4)

plotvect = CD8pos_tumor(:,1) + (rand(size(CD8pos_tumor(:,1)))*0.4) - 0.2;
mean_forbar = [];
for i = 1:2
    mean_forbar(i) = nanmean(CD8pos_tumor(CD8pos_tumor(:,1)==i,4));
end
bar([1:size(forplot.group_mat,1)]+[0.375 -0.375],mean_forbar,1,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
hold on
scatter(plotvect,CD8pos_tumor(:,4),4,CD8pos_tumor(:,1),'filled')
caxis([0.5 size(forplot.group_mat,1)])
xlim([0.5 size(forplot.group_mat,1)+.5])
ylim([0 max(CD8pos_tumor(:,4))*1.1])
colormap('jet')
title(['Frac GzmB+ CD8 cells per tumor'])

%%%% print csv's
nameforfile = ['20201130_KP_Vaccine_Mouse_GzmBPosFrac_cytoTcells_alllung.csv'];
forprint = sortrows([forplot.group_mat(:) CD8pos_total_mouse(:)]);
csvwrite(nameforfile,forprint)

nameforfile = ['20201130_KP_Vaccine_Mouse_GzmBPosFrac_cytoTcells_intumor.csv'];
forprint = sortrows([forplot.group_mat(:) CD8pos_tumor_mouse(:)]);
csvwrite(nameforfile,forprint)

nameforfile = ['20201130_KP_Vaccine_Mouse_GzmBPosFrac_cytoTcells_pertumor.csv'];
forprint = CD8pos_tumor;
csvwrite(nameforfile,forprint)

%% Ki67 in CD8+ exploration
probe_round = find(strcmp(options.Markers,'Ki67'));

% find cd8 cells in the whole mouse and in the tumors
CD8cells = CellType.Matrix(:,4) == 1113 & Filter.all(:,probe_round)>0 & ROIResults.SpleenIndex == 0;
CD8cells_intumor = CellType.Matrix(:,4) == 1113 & ROIResults.TumorIndex > 0 & Filter.all(:,probe_round)>0;

linecols = {'b','r'};
clear temp

CD8pos_total_mouse = zeros(size(forplot.group_mat)) + NaN;
CD8pos_tumor_mouse = zeros(size(forplot.group_mat)) + NaN;
CD8pos_tumor = [];

count = 0;

for g = 1:size(forplot.group_mat,1)
    gzmB_cum = [];
    tum_ind_cum = [];
    
    % loop around the mouse within each condition
    for i = 1:size(forplot.group_mat,2)
        % isolate cells of specific mouse
        rois_ind = find(options.MouseNum==forplot.group_mat(g,i));
        clear D
        temp_CD8cells = 0;
        if length(rois_ind) == 0
            gzmb_pos(g,i) = NaN;
            continue
        end
        
        cd8_mouse = sum(CD8cells & MorpResults.Indexes == rois_ind);
        cd8pos_mouse = sum(CD8cells & MorpResults.Indexes == rois_ind & NormResults.MedianNucNorm(:,probe_round)>0);
        
        cd8_tumor = sum(CD8cells_intumor & MorpResults.Indexes == rois_ind);
        cd8pos_tumor = sum(CD8cells_intumor & MorpResults.Indexes == rois_ind & NormResults.MedianNucNorm(:,probe_round)>0);
        
        CD8pos_total_mouse(g,i) = cd8pos_mouse/cd8_mouse ; 
        CD8pos_tumor_mouse(g,i) = cd8pos_tumor/cd8_tumor ;  
        
        
        % get stats per tumor
        tumor_indices = unique(ROIResults.TumorIndex(MorpResults.Indexes == rois_ind));
        
        for t = 2:length(tumor_indices)
            % find tumor cells - skip the first index as it is == 0 ie not
            % tumor lung tissue
            tum_ind = ROIResults.TumorIndex == tumor_indices(t) & MorpResults.Indexes == rois_ind ;
            % calculate num of CD8 and CD* positive for the proeb
            CD8_pertum = sum(CD8cells & tum_ind);
            CD8pos_pertum = sum(CD8cells & tum_ind & NormResults.MedianNucNorm(:,probe_round)>0);
            
            % check that there are CD8 in this tumor
            if CD8_pertum == 0
                continue
            end
            
            count = count + 1;
            CD8pos_tumor(count,1:4) = [g forplot.group_mat(g,i) tumor_indices(t) CD8pos_pertum/CD8_pertum] ;
        end
        
        gzmb_CD8cell = NormResults.MedianNucNorm(CD8cells & MorpResults.Indexes == rois_ind,probe_round);
        tumind_temp = ROIResults.TumorIndex(CD8cells & MorpResults.Indexes == rois_ind) > 0;
        
        gzmB_cum = [gzmB_cum; gzmb_CD8cell];
        tum_ind_cum = [tum_ind_cum ; tumind_temp];
        
    end
    
    figure(34)
    subplot(2,2,1)
    [n,h]=ksdensity(gzmB_cum);
    plot(h,n,linecols{g},'LineWidth',2)
    hold on
    xlim([-2500 2500])
    
    GrzBpos_prop = sum(gzmB_cum > 0)/length(gzmB_cum);
    
    if g == 1
        text(1000,(1)*10^(-3),'Overall Ki67+ %' ,'HorizontalAlignment','center')
    end
    text(1000,(1-g*0.05)*10^(-3),[forplot.grouplist{g} ' =  ' num2str(round(GrzBpos_prop*100,2)) '%'],'HorizontalAlignment','center')
    
    title(['Single cell Ki67 distribution by group'])
    xlabel('Normalized Ki67 CyCIF signal')
    ylabel('ks pdf (frequency)')
    
    GrzBpos_prop_loc(g,1) = sum(gzmB_cum > 0)/length(gzmB_cum);
    GrzBpos_prop_loc(g,2) = sum(gzmB_cum(tum_ind_cum > 0) > 0)/length(gzmB_cum(tum_ind_cum > 0));
    
    
end

figure(34)
subplot(2,2,1)
legend(forplot.grouplist)

subplot(2,2,2)
bar(GrzBpos_prop_loc')
legend(forplot.grouplist)
set(gca,'XTick',1:size(forplot.group_mat,1),'XTickLabel',{'ALL','IN'})%,'XTickLabelRotation',45)

subplot(2,2,3)
mouse_stat_plot([CD8pos_total_mouse; CD8pos_tumor_mouse])
labels = {['All ' forplot.grouplist{1}],['All ' forplot.grouplist{2}],['IN ' forplot.grouplist{1}],['IN ' forplot.grouplist{2}]};
set(gca,'XTick',1:size(forplot.group_mat,1)*2,'XTickLabel',labels)
[~,p_All]=ttest2(CD8pos_total_mouse(1,:),CD8pos_total_mouse(2,:));
[~,p_IN]=ttest2(CD8pos_tumor_mouse(1,:),CD8pos_total_mouse(2,:));
title(['p-All = ' num2str(round(p_All,5)) ' & p-IN = ' num2str(round(p_IN,5))])

subplot(2,2,4)

plotvect = CD8pos_tumor(:,1) + (rand(size(CD8pos_tumor(:,1)))*0.4) - 0.2;
mean_forbar = [];
for i = 1:2
    mean_forbar(i) = nanmean(CD8pos_tumor(CD8pos_tumor(:,1)==i,4));
end
bar([1:size(forplot.group_mat,1)]+[0.375 -0.375],mean_forbar,1,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
hold on
scatter(plotvect,CD8pos_tumor(:,4),4,CD8pos_tumor(:,1),'filled')
caxis([0.5 size(forplot.group_mat,1)])
xlim([0.5 size(forplot.group_mat,1)+.5])
ylim([0 max(CD8pos_tumor(:,4))*1.1])
colormap('jet')
title(['Frac Ki67+ CD8 cells per tumor'])

%%%% print csv's
nameforfile = ['20201130_KP_Vaccine_Mouse_Ki67PosFrac_cytoTcells_alllung.csv'];
forprint = sortrows([forplot.group_mat(:) CD8pos_total_mouse(:)]);
csvwrite(nameforfile,forprint)

nameforfile = ['20201130_KP_Vaccine_Mouse_Ki67PosFrac_cytoTcells_intumor.csv'];
forprint = sortrows([forplot.group_mat(:) CD8pos_tumor_mouse(:)]);
csvwrite(nameforfile,forprint)

nameforfile = ['20201130_KP_Vaccine_Mouse_Ki67PosFrac_cytoTcells_pertumor.csv'];
forprint = CD8pos_tumor;
csvwrite(nameforfile,forprint)

%% Tumor burden quantification
scale = 8;
plotcheck = 0;
cum_mat = [];
count  = 0;
tum_count = zeros(size(forplot.group_mat))+NaN;
tum_area = tum_count;
tum_area_corrected = tum_count;

for g = 1:size(forplot.group_mat,1)
    tum_index = [];
    
    % loop around the mouse within each condition
    for i = 1:size(forplot.group_mat,2)
        rois_ind = find(options.MouseNum==forplot.group_mat(g,i));
        if length(rois_ind) == 0
            continue
        end
        % isolate cells of specific mouse
        rois_ind = find(options.MouseNum==forplot.group_mat(g,i));
        allcells = MorpResults.Indexes == rois_ind & ROIResults.SpleenIndex == 0;
        tumcells = MorpResults.Indexes == rois_ind & ROIResults.TumorIndex > 0;
        
        X = MorpResults.X(allcells);        
        Y = MorpResults.Y(allcells);
        
        X_tum = MorpResults.X(tumcells);        
        Y_tum = MorpResults.Y(tumcells);
        tum_num = ROIResults.TumorIndex(tumcells);
        
        lung_img = uint8(zeros(ceil(max(X)/scale),ceil(max(Y)/scale)));
        tum_img = uint8(zeros(ceil(max(X)/scale),ceil(max(Y)/scale)));
        
        cell_ind = sub2ind(size(lung_img),round(X/scale),round(Y/scale));
        tum_ind = sub2ind(size(tum_img),round(X_tum/scale),round(Y_tum/scale));
        
        lung_img(cell_ind) = 1;
        lung_img = imdilate(lung_img,strel('octagon',15));
        lung_img = imfill(lung_img,'holes');
        lung_img = imerode(lung_img,strel('octagon',15));
        
        tum_img(tum_ind) = tum_num;
        tum_img = imdilate(tum_img,strel('octagon',6));
        tum_img = imfill(tum_img,'holes');
        tum_img = imerode(tum_img,strel('octagon',6));
        tum_img(lung_img==0) = 0;
        
        if plotcheck
            figure
            subplot(1,2,1)
            imshow(lung_img,[])
            subplot(1,2,2)
            imshow(tum_img,[])
        end
        
        % calculate area of lung and tumor
        LungStats = regionprops(lung_img,'Area');
        TumStats = regionprops(tum_img,'Area');

        for tum = 1:max(ROIResults.TumorIndex(tumcells))
            
            tum_index = tumcells & ROIResults.TumorIndex == tum;
            
            if sum(tum_index) == 0
                continue
            end
            
            count = count + 1;
            totcount = sum(tum_index);
            epicount = sum(CellType.Matrix(tum_index,1)==2);
            immcount = sum(CellType.Matrix(tum_index,1)==1);

            cum_mat(count,1:7) = [g forplot.group_mat(g,i) totcount epicount TumStats(tum).Area TumStats(tum).Area*(epicount/totcount) LungStats.Area];
        end
 
    end
end%%

for g = 1:size(forplot.group_mat,1)
    % loop around the mouse within each condition
    for i = 1:size(forplot.group_mat,2)
        if forplot.group_mat(g,i) == 0
            continue
        end
        sub_ind = cum_mat(:,1)==g & cum_mat(:,2)==forplot.group_mat(g,i);
        tum_area(g,i) =           mean(cum_mat(sub_ind,5)) / mean(cum_mat(sub_ind,7)); 
        tum_area_corrected(g,i) = mean(cum_mat(sub_ind,6)) / mean(cum_mat(sub_ind,7)); 
        tum_count(g,i) = sum(sub_ind);
    end
end

figure(61)

subplot(1,3,1)
mouse_stat_plot(tum_count)
[~,p]=ttest2(tum_count(1,:),tum_count(2,:));
title(['p=' num2str(round(p,5))])
set(gca,'XTick',1:size(forplot.group_mat,1),'XTickLabel',forplot.grouplist)
ylabel('Total Number of Tumors')

subplot(1,3,2)
mouse_stat_plot(tum_area)
[~,p]=ttest2(tum_area(1,:),tum_area(2,:));
title(['p=' num2str(round(p,5))])
set(gca,'XTick',1:size(forplot.group_mat,1),'XTickLabel',forplot.grouplist)
ylabel('Mean Tumor Area')

subplot(1,3,3)
mouse_stat_plot(tum_area_corrected)
[~,p]=ttest2(tum_area_corrected(1,:),tum_area_corrected(2,:));
title(['p=' num2str(round(p,5))])
set(gca,'XTick',1:size(forplot.group_mat,1),'XTickLabel',forplot.grouplist)
ylabel('Mean Tumor Corrected')

suptitle('Tumor Burden Calculations - based on area')

nameforfile = ['20201130_KP_Vaccine_Mouse_Tumor_Burden_byArea.csv'];
forprint = cum_mat;
csvwrite(nameforfile,forprint)

%% Tumor burden quantification by cell count
tot_cells = [];
tum_cells = [];
epi_cells = [];
imm_cells = [];

cum_mat = [];
mean_tum_count = zeros(2,8) + NaN;
mean_epi_count = mean_tum_count;
tum_count = mean_tum_count;
count = 0;

for g = 1:size(forplot.group_mat,1)
    tum_index = [];
    
    % loop around the mouse within each condition
    for i = 1:size(forplot.group_mat,2)
        % isolate cells of specific mouse
        rois_ind = find(options.MouseNum==forplot.group_mat(g,i));
        clear D
        temp_CD8cells = 0;
        if length(rois_ind) == 0
            tot_cells(g,i) = NaN;
            tum_cells(g,i) = NaN;
            epi_cells(g,i) = NaN;
            imm_cells(g,i) = NaN;
            continue
        end
        for r = 1:length(rois_ind)
            index = MorpResults.Indexes == rois_ind(r) & ROIResults.SpleenIndex == 0;
            
            tot_cells(g,i) = sum(index);
            tum_cells(g,i) = sum(index & ROIResults.TumorIndex > 0);
            epi_cells(g,i) = sum(index & ROIResults.TumorIndex > 0 & CellType.Matrix(:,1) == 2);
            imm_cells(g,i) = sum(index & ROIResults.TumorIndex > 0 & CellType.Matrix(:,1) == 1);
            
            for tum = 1:max(ROIResults.TumorIndex(index))
                count = count + 1;
                tum_index = ismember(MorpResults.Indexes,rois_ind) & ROIResults.SpleenIndex == 0 & ROIResults.TumorIndex == tum;

                totcount = sum(tum_index);
                epicount = sum(CellType.Matrix(tum_index,1)==2);
                
                cum_mat(count,1:5) = [g forplot.group_mat(g,i) totcount epicount tot_cells(g,i)];
            end
            
        end
    end
end

for g = 1:size(forplot.group_mat,1)
    % loop around the mouse within each condition
    for i = 1:size(forplot.group_mat,2)
        mouse = forplot.group_mat(g,i);
        if mouse == 0
            continue
        end
        mean_tum_count(g,i) = mean(cum_mat(cum_mat(:,1)==g & cum_mat(:,2)==mouse,3)) / tot_cells(g,i); 
        mean_epi_count(g,i) = mean(cum_mat(cum_mat(:,1)==g & cum_mat(:,2)==mouse,4)) / tot_cells(g,i); 
        tum_count(g,i) = sum(cum_mat(:,1)==g & cum_mat(:,2)==mouse); 
    end
end

figure(71)

subplot(1,3,1)
mouse_stat_plot(tum_count)
[~,p]=ttest2(tum_count(1,:),tum_count(2,:));
title(['p=' num2str(round(p,5))])
set(gca,'XTick',1:size(forplot.group_mat,1),'XTickLabel',forplot.grouplist)
ylabel('Total Number of Tumors')

subplot(1,3,2)
mouse_stat_plot(mean_tum_count)
[~,p]=ttest2(mean_tum_count(1,:),mean_tum_count(2,:));
title(['p=' num2str(round(p,5))])
set(gca,'XTick',1:size(forplot.group_mat,1),'XTickLabel',forplot.grouplist)
ylabel('Mean Tumor CellNum / Total Cells')

subplot(1,3,3)
mouse_stat_plot(mean_epi_count)
[~,p]=ttest2(mean_epi_count(1,:),mean_epi_count(2,:));
title(['p=' num2str(round(p,5))])
set(gca,'XTick',1:size(forplot.group_mat,1),'XTickLabel',forplot.grouplist)
ylabel('Mean Epith Tumor CellNum / Total Cells')

suptitle('Tumor Burden Calculations - based on cellcount')

nameforfile = ['20201130_KP_Vaccine_Mouse_Tumor_Burden_byCellcount.csv'];
forprint = cum_mat;
csvwrite(nameforfile,forprint)

%% B cells explorationa
ki67_round = find(strcmp(options.Markers,'Ki67'));
bcells = CellType.Matrix(:,3) == 112 & Filter.all(:,ki67_round)>0;
linecols = {'b','r'};
clear temp
d = 8;

d_thresh = log2(100); % cut off the distance to 100 pixels

for g = 1:size(forplot.group_mat,1)
    % loop around the mouse within each condition
    ki_cum = [];
    D_cum = [];
    for i = 1:size(forplot.group_mat,2)
        % isolate cells of specific mouse
        rois_ind = find(options.MouseNum==forplot.group_mat(g,i));
        clear D
        temp_bcells = 0;
        if length(rois_ind) == 0
            continue
        end
        for r = 1:length(rois_ind)
            index = MorpResults.Indexes == rois_ind(r) & ROIResults.SpleenIndex == 0;
            
            % count B cells
            temp_bcells = temp_bcells + sum(index & bcells);
            
            cent_x = double(ceil(MorpResults.X(index & bcells)));
            cent_y = double(ceil(MorpResults.Y(index & bcells)));
            ki67_bcell = NormResults.MedianNucNorm(index & bcells,find(strcmp(options.Markers,'Ki67')));
            ki67pos(g,i) = sum(ki67_bcell>0)/length(ki67_bcell);
            
            [~,D] = knnsearch([cent_x cent_y],[cent_x cent_y],'K',d);
        end
        num_bcells(g,i) = temp_bcells;
      
        figure(50+g)
        subplot(2,2,2)
        [n,h]=ksdensity(log2(D(:,d)+1));
        plot(h,n,'color',[0.2 0.5 1])
        hold on
        
        subplot(2,2,3)
        [n,h]=ksdensity(ki67_bcell);
        plot(n,h,'color',[0.2 0.5 1])
        hold on

        ki_cum = [ki_cum; ki67_bcell];
        D_cum  = [D_cum; log2(D(:,d)+1)];
        
    end
    figure(50+g)
    subplot(2,2,4)
    scatter(D_cum,ki_cum,10,'b','filled','MarkerFaceAlpha',0.05)
    axis([ 4.5 10.5 -1000 2000 ])
    hold on
    
    subplot(2,2,2)
    [n,h]=ksdensity(D_cum);
    plot(h,n,linecols{g},'LineWidth',3)
    xlim([4.5 10.5])
    title(['log2 Distance from ' num2str(d) '^t^h B cells'])
    
    subplot(2,2,3)
    [n,h]=ksdensity(ki_cum);
    plot(n,h,linecols{g},'LineWidth',3)
    ylim([-1000 2000])
    set(gca,'XDir','reverse')
    title('Ki67')
    
    suptitle(['B Cells ' forplot.grouplist{g}])
    
    figure(53)
    [n,h]=ksdensity(ki_cum(D_cum>d_thresh));
    plot(h,n,linecols{g},'LineWidth',2)
    hold on
    [n,h]=ksdensity(ki_cum(D_cum<=d_thresh));
    plot(h,n,linecols{g},'LineWidth',2,'LineStyle','--')
    xlim([-1000 2000])
    title(['Ki67 binned by distance from ' num2str(d) '^t^h B cells with ' num2str(d_thresh) ' pix thresh'])
    
    
end

figure(53)
legend({[forplot.grouplist{1} ' D > th'],[forplot.grouplist{1} ' D < th'],[forplot.grouplist{2} ' D > th'],[forplot.grouplist{2} ' D < th']})

%%
% close all
axis_max = zeros(1,length(CellType.Classes));
norm_cellnum = cell(2,size(forplot.group_mat,1));

% loop around either the whole lung or the tumor regions only (as a binary)
for loc = 1:2
    % loop around the mouse cohort/conditions
    for g = 1:size(forplot.group_mat,1)
        % loop around the mouse within each condition
        for i = 1:size(forplot.group_mat,2)
            disp([num2str(g) ' ' num2str(i)])
            rois_ind = find(options.MouseNum==forplot.group_mat(g,i));
            index = ismember(MorpResults.Indexes,rois_ind) & ROIResults.SpleenIndex == 0 & ROIResults.TumorIndex > loc-2;

            for type = 1:max(CellType.index)
                % CHOOSE NORMALIZATION to either total cells
                norm_cellnum{loc,type}(g,i) = sum(CellType.Matrix(index,CellType.layer(type))==CellType.codes(type))/sum(CellType.Matrix(index,1)>0)*100;
                % or to only immmune cells
%                 norm_cellnum{loc,type}(g,i) = sum(CellType.Matrix(index,CellType.layer(type))==CellType.codes(type))/sum(CellType.Matrix(index,1)==1)*100;
                % or to lymphocyte cells
%                 norm_cellnum{loc,type}(g,i) = sum(CellType.Matrix(index,CellType.layer(type))==CellType.codes(type))/sum(CellType.Matrix(index,2)==11)*100;

                if CellType.codes(type) == 1113 % ie cytotoxic t cells
                    ind_Tcytox = CellType.Matrix(index,CellType.layer(type))==CellType.codes(type);
                    labels_Tcytox = {'T cytox','GzmB','Perforin','Ki67','PD-1','Tim3','PD1+Tim3+'};
                    ind_sub{1} = ones(size(ind_Tcytox));
                    ind_sub{2} = NormResults.MedianNucNorm(index,find(strcmp(labels_Tcytox{2},options.Markers))) > 0;
                    ind_sub{3} = NormResults.MedianNucNorm(index,find(strcmp(labels_Tcytox{3},options.Markers))) > 0;
                    ind_sub{4} = NormResults.MedianNucNorm(index,find(strcmp(labels_Tcytox{4},options.Markers))) > 0;
                    ind_sub{5} = NormResults.MedianCytNorm(index,find(strcmp(labels_Tcytox{5},options.Markers))) > 0;
                    ind_sub{6} = NormResults.MedianCytNorm(index,find(strcmp(labels_Tcytox{6},options.Markers))) > 0;
                    ind_sub{7} = NormResults.MedianCytNorm(index,find(strcmp(labels_Tcytox{5},options.Markers))) > 0 & ...
                                 NormResults.MedianCytNorm(index,find(strcmp(labels_Tcytox{6},options.Markers))) > 0;

                    for sub_i = 1:length(ind_sub)
                        % CHOOSE NORMALIZATION to either total cells
%                         norm_subcellnum{loc,sub_i}(g,i) = sum(ind_Tcytox & ind_sub{sub_i})/sum(CellType.Matrix(index,1)>0)*100;
                        % or to only immmune cells
%                         norm_subcellnum{loc,sub_i}(g,i) = sum(ind_Tcytox & ind_sub{sub_i})/sum(CellType.Matrix(index,1)==1)*100;
                        % or to only immmune cells
                        norm_subcellnum{loc,sub_i}(g,i) = sum(ind_Tcytox & ind_sub{sub_i})/sum(ind_Tcytox)*100;
                    end
                end

            end
        end
    end
end

plotvect = [ones(5,1)-(6-(1:5)'-3)*0.1; 2*ones(5,1)-(6-(1:5)'-3)*0.1; 3*ones(5,1)-(6-(1:5)'-3)*0.1; 4*ones(5,1)-(6-(1:5)'-3)*0.1];

for loc = 1:2
    figure(loc+100)
    for type = 1:max(CellType.index)
        temp_mat = norm_cellnum{loc,type}';
        median_forbar = median(temp_mat);
        subplot(4,4,type)
        bar([1:4],median_forbar,0.5,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
        hold on
        scatter(plotvect,temp_mat(:),20,round(plotvect),'filled')
        caxis([0 4.5])
        xlim([0.5 4.5])
        ylim([0 max(temp_mat(:))*1.1])
        colormap('jet')
        set(gca,'XTick',1:4,'XTickLabel',forplot.grouplist,'XTickLabelRotation',45)
        title(CellType.names{type})
        if mod(type-1,4)==0 && loc == 1
            ylabel('% of total cells')
        elseif mod(type-1,4)==0 && loc == 2
            ylabel('% of tumor cells')
        end
    end

    figure(loc+200)
    for type = 1:length(ind_sub)
        temp_mat = norm_subcellnum{loc,type}';
        median_forbar = median(temp_mat);
        subplot(2,4,type)
        bar([1:4],median_forbar,0.5,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
        hold on
        scatter(plotvect,temp_mat(:),20,round(plotvect),'filled')
        caxis([0 4.5])
        xlim([0.5 4.5])
        ylim([0 max(temp_mat(:))*1.1])
        colormap('jet')
        set(gca,'XTick',1:4,'XTickLabel',forplot.grouplist,'XTickLabelRotation',45)
        title(labels_Tcytox{type})
        if mod(type-1,4)==0 && loc == 1
            ylabel('% of total cells')
        elseif mod(type-1,4)==0 && loc == 2
            ylabel('% of tumor cells')
        end
    end
end 

figure(101)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.3 0.8 ])
pause(3)
suptitle('Whole Lung Area')

figure(102)
set(gcf,'Units','normalized','Position',[0.1 0.1 0.32 0.8 ])
pause(3)
suptitle('Tumor ROIs Area')

figure(201)
suptitle('Whole Lung Area')

figure(202)
suptitle('Tumor ROIs Area')
