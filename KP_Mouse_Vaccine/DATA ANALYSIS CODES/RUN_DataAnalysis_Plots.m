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
load([analfolder resufolder 'Results_Nets_' date '_dist25.mat'])
options.figOpt = 0;

filename.basefolder = basefolder;
filename.analfolder = analfolder;
filename.resufolder = resufolder; 
options.date = date;

save([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'],'filename','options')

forplot.grouplist = {'ctrl','vax'};
forplot.group_mat = [1117:1123 0;  1124:1131];%
% forplot.group_mat = [1117:1123;   1126:1131 0;];%
% forplot.group_mat = [1124 1126 1127 1130;  1125 1128 1129 1131];
% forplot.group_mat = [1117:1123; 1124 1126 1127 1130 0 0 0];
% forplot.group_mat = [1117:1123; 1125 1128 1129 1131 0 0 0];

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

%%% compare cell density to lymphonets networks
for i = 1:length(CellType.names)
    figure(100)
    subplot(4,4,i)
%     plot(nets_num([test_group1 test_group2],:)',D_mean{i}([test_group1 test_group2],:)','.','MarkerSize',25)
    plot(nets_num([test_group1 ],:)',D_mean{i}([test_group1 ],:)','.','MarkerSize',25)
    hold on
    title(CellType.names{i})
    if i == 1
        legend(forplot.grouplist)
    end
%     x = nets_num([test_group1 test_group2],:);
%     x = x(:);
%     y = D_mean{i}([test_group1 test_group2],:);
%     y = y(:);
    
    x = nets_num([test_group1 ],:);
    x = x(:);
    y = D_mean{i}([test_group1 ],:);
    y = y(:);
    
    x = x(~isnan(y));
    y = y(~isnan(y));
    
    [p,S]=polyfit(x,y,1);
    R2 = 1 - (S.normr/norm(y - mean(y)))^2;
    plot([min(x); max(x)],p(2)+p(1)*([min(x); max(x)]),'--k')
    title([CellType.names{i} ' - R^2=' num2str(round(R2,4))])
    xlabel(['# nets > ' num2str(th_nets_cellnum) ' cells'])
    ylabel('Mean Cell Density')
    axis([0 max(x)*1.2+0.0001 0 max(y)*1.2+0.0001])
end

set(gcf,'Units','normalized','Position',[0.1 0.1 0.5 0.8 ])
pause(3)
suptitle('Cell Retention by Lymphonets')
         

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
%     [n,h]=ksdensity(gzmB_cum(tum_ind_cum == 0));
%     plot(h,n,linecols{g},'LineWidth',2)
%     hold on
%     [n,h]=ksdensity(gzmB_cum(tum_ind_cum > 0));
%     plot(h,n,linecols{g},'LineWidth',2,'LineStyle','--')
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
%     [n,h]=ksdensity(gzmB_cum(tum_ind_cum == 0));
%     plot(h,n,linecols{g},'LineWidth',2)
%     hold on
%     [n,h]=ksdensity(gzmB_cum(tum_ind_cum > 0));
%     plot(h,n,linecols{g},'LineWidth',2,'LineStyle','--')
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


%% calculate summary metric for the nets

allnets.size  = [];
allnets.bfrac = [];
allnets.tfrac = [];
allnets.tRfrac = [];
allnets.tHfrac = [];
allnets.tCfrac = [];
allnets.mouse = [];
allnets.mousegroup = [];
for j = 1:length(options.MouseNum)
    allnets.size = [allnets.size; LymphoNets.Summary{j}.NumCells'];
    allnets.bfrac = [allnets.bfrac; LymphoNets.Summary{j}.Bcells'./LymphoNets.Summary{j}.NumCells'];
    allnets.tfrac = [allnets.tfrac; LymphoNets.Summary{j}.Tcells'./LymphoNets.Summary{j}.NumCells'];
    allnets.tRfrac = [allnets.tRfrac; LymphoNets.Summary{j}.Tregs'./LymphoNets.Summary{j}.NumCells'];
    allnets.tHfrac = [allnets.tHfrac; LymphoNets.Summary{j}.Thelps'./LymphoNets.Summary{j}.NumCells'];
    allnets.tCfrac = [allnets.tCfrac; LymphoNets.Summary{j}.Tcytox'./LymphoNets.Summary{j}.NumCells'];
    allnets.mouse = [allnets.mouse; zeros(length(LymphoNets.Summary{j}.NumCells),1)+options.MouseNum(j)];
    allnets.mousegroup = [allnets.mousegroup; zeros(length(LymphoNets.Summary{j}.NumCells),1)+options.MouseGroup(j)];
    
end

colorbygroup = {'b','c','r','m','g'};
cum_nets = cell(1,2);
for j = 1:length(options.MouseNum)
    index = MorpResults.Indexes == j;
    if options.MouseGroup(j) ==  5 % < 5
        
        figure(5000+j)
        scatter(MorpResults.X(index & ROIResults.TumorIndex > 0),MorpResults.Y(index  & ROIResults.TumorIndex > 0 ),10,[255 153 51]/255,'filled')
        set(gca,'Ydir','reverse')
        hold on
        scatter(MorpResults.X(index & ROIResults.TumorIndex == 0),MorpResults.Y(index  & ROIResults.TumorIndex == 0 ),10,[51 255 255]/255,'filled','MarkerFaceAlpha',1)%0.05)
        scatter(MorpResults.X(index & LymphoNets.NetworkID > 0),MorpResults.Y(index  & LymphoNets.NetworkID > 0),10,[255 51 255]/255,'filled')
        
        viscircles([LymphoNets.Summary{j}.X_mean' LymphoNets.Summary{j}.Y_mean'],LymphoNets.Summary{j}.NumCells'+0.1,'Color','k','LineStyle','--');
        viscircles([LymphoNets.Summary{j}.X_mean' LymphoNets.Summary{j}.Y_mean'],LymphoNets.Summary{j}.Tcells'+0.1,'Color','b');
        viscircles([LymphoNets.Summary{j}.X_mean' LymphoNets.Summary{j}.Y_mean'],LymphoNets.Summary{j}.Bcells'+0.1,'Color','r');
       
    end
    cum_nets{options.MouseGroup(j)} = [cum_nets{options.MouseGroup(j)} LymphoNets.Summary{j}.NumCells];
end


th_cell = 5;

for i = 1:2
    mean(cum_nets{i}(cum_nets{i}>th_cell))
    sum(cum_nets{i}>th_cell)
    figure(1001)
    subplot(1,2,i)
    [n,h]=hist(cum_nets{i}(cum_nets{i}>th_cell),linspace(th_cell,th_cell*12,11));
    bar(h,n)%/sum(n))
    hold on
%     ylim([0 max([80 n])*1.02])
    ylim([0 300])
    title(forplot.grouplist{i})
    xlabel('Network size')
end
suptitle('Distribution of nets by size')

%%% COMPOSITION
mean_Bfrac = [];
mean_Tfrac = [];
mean_tCfrac = [];
mean_tHfrac = [];
mean_tRfrac = [];
tot= [];
x = [];
ylim_cum = 0.7;

bins = 40;
jump = 5;
max_net = bins*jump;
% find
for j = 1:2
    for i = 1:bins
        % find nets of certain size
        ind_size = allnets.size > jump*(i-1) & allnets.size <= jump*i & allnets.mousegroup == j;
        mean_Bfrac(i,:) = [nanmean(allnets.bfrac(ind_size)) prctile(allnets.bfrac(ind_size),25) prctile(allnets.bfrac(ind_size),75) ]; 
        mean_Tfrac(i,:) = [nanmean(allnets.tfrac(ind_size)) prctile(allnets.tfrac(ind_size),25) prctile(allnets.tfrac(ind_size),75) ]; 
        mean_tCfrac(i,:) = [nanmean(allnets.tCfrac(ind_size)./allnets.tfrac(ind_size)) prctile(allnets.tCfrac(ind_size)./allnets.tfrac(ind_size),25) prctile(allnets.tCfrac(ind_size)./allnets.tfrac(ind_size),75) ]; 
        mean_tHfrac(i,:) = [nanmean(allnets.tHfrac(ind_size)./allnets.tfrac(ind_size)) prctile(allnets.tHfrac(ind_size)./allnets.tfrac(ind_size),25) prctile(allnets.tHfrac(ind_size)./allnets.tfrac(ind_size),75) ]; 
        mean_tRfrac(i,:) = [nanmean(allnets.tRfrac(ind_size)./allnets.tfrac(ind_size)) prctile(allnets.tRfrac(ind_size)./allnets.tfrac(ind_size),25) prctile(allnets.tRfrac(ind_size)./allnets.tfrac(ind_size),75) ];  
        x(i,1) = i*jump - jump/2;
        tot(i) = nansum(ind_size);
        if tot(i) < 10
            mean_Bfrac(i,:) = [ NaN NaN NaN];
            mean_Tfrac(i,:) = [ NaN NaN NaN];
            mean_tCfrac(i,:) = [ NaN NaN NaN];
            mean_tHfrac(i,:) = [ NaN NaN NaN];
            mean_tRfrac(i,:) = [ NaN NaN NaN];
        end
    end
    figure(1001+j)
    subplot(1,2,1)
    scatter(x-0.5,mean_Bfrac(:,1),log2(tot+1.1)*3,'r','filled')
    hold on
    scatter(x+0.5,mean_Tfrac(:,1),log2(tot+1.1)*3,'b','filled')
    plot([x-.5 x-.5]',mean_Bfrac(:,2:3)','-r')
    
    plot([x+.5 x+.5]',mean_Tfrac(:,2:3)','-b')
    legend('B cell frac','T cell frac')
    ylim_cum = max([ylim_cum max([mean_Bfrac,mean_Tfrac])]);
    ylim([0 ylim_cum*1.05])
    xlabel('Network size')
    ylabel('Fraction of network')
    
    subplot(1,2,2)
    scatter(x,mean_tCfrac(:,1),5*3,'filled')
    hold on
    scatter(x,mean_tHfrac(:,1),5*3,'filled')
    scatter(x,mean_tRfrac(:,1),5*3,'filled')
    legend('Tc cell frac','Th cell frac','Tr cell frac')
    ylim([0 ylim_cum*1.05])
    xlabel('Network size')
    ylabel('Fraction within T cells')
    
    suptitle(forplot.grouplist{j})
end

% %% does the location of T cell phenotypic markers change in or out of network?
% 
% marker_list = {'GranzB','PD-1','Tim3-2','Ki67','Perforin'};
% thresh =      [ 200      0      0        0      200];%[ 300      200    200      200];
% celltype_list = {'B','T cytotox','T reg','T helper'};%,'T reg','T helper'};
% treat_list = {'Luc','LucOS'};
% 
% th_cell = 15;
% for t = 1:length(celltype_list)
%     celltype_code = CellType.codes(find(strcmp(CellType.names,celltype_list{t})));
%     celltype_layer = ceil(log10(celltype_code+1));
%     for g = 1:5
%         if g < 5
%             rois_ind = find(ismember(options.MouseNum,forplot.group_mat(g,:)));
%             index = ROIResults.SpleenIndex == 0 & CellType.Matrix(:,celltype_layer) == celltype_code & ismember(MorpResults.Indexes,rois_ind);
%         else
%             index = ROIResults.SpleenIndex == 0 & CellType.Matrix(:,celltype_layer) == celltype_code;
%         end
%         index_in = index & LymphoNets.Size > th_cell;
%         index_out =  index & LymphoNets.Size <= th_cell;
% 
%         for m = 1:length(marker_list)
%             mark = find(strcmp(options.Markers,marker_list{m}));
% 
% 
%             vect = NormResults.MedianNucNorm(index,mark);
%             vect_in = NormResults.MedianNucNorm(index_in,mark);
%             vect_out = NormResults.MedianNucNorm(index_out,mark);
% 
%             figure(3000+t*10+g)
%             if strcmp(marker_list{m},'PD-1')
%                 med_frac_in = sum(vect_in > thresh(m)) / length(vect_in);
%                 med_frac_out = sum(vect_out > thresh(m)) / length(vect_out);
%                 high_frac_in = sum(vect_in > thresh(m)+750) / length(vect_in);
%                 high_frac_out = sum(vect_out > thresh(m)+750) / length(vect_out);
%                 subplot(2,length(marker_list),length(marker_list)+m)
%                 bar([med_frac_out med_frac_in; high_frac_out high_frac_in])
%                 xlabel([num2str(length(vect_in))])
%             else
%                 pos_frac_in = sum(vect_in > thresh(m)) / length(vect_in);
%                 pos_frac_out = sum(vect_out > thresh(m)) / length(vect_out);
%                 subplot(2,length(marker_list),length(marker_list)+m)
%                 bar([pos_frac_out pos_frac_in])
%                 xlabel([num2str(length(vect_in))])
%             end
%             figure(3000+t*10+g)
%             subplot(2,length(marker_list),m)
%             ksdensity(vect_in)
%             hold on
%             ksdensity(vect_out)
%             xlim([-1500 1500])
%             title(marker_list{m})
% 
% 
% 
%         end 
%         figure(3000+t*10+g)
%         try
%             suptitle([celltype_list{t} ' - ' forplot.grouplist{g}])
%         catch
%             suptitle([celltype_list{t} ' - All mouse groups'])
%         end
%     end
% end
% 
% %%
% 
% clear totcells epicells prolifer immune tcells
% 
% axis_max = zeros(1,length(CellType.Classes));
% grouplist =  forplot.grouplist;
% group_mat = forplot.group_mat;
% 
% % norm_cellnum = cell(size(group_mat));
% 
% dotsizemult = 3;
% maxmouse = 20;
% colorlist = {'b','r'};
% label_all = {};
% mat = [];
% nontumor = [];
% label_mat = {'Mouse','Totcells','Epithcells','T','B','NK','Macro','Dendr','Neutro','T reg','T helper','T cytox','T cytox GranzB+',...
%              'T cytox Perforin+','T cytox Ki67+','T cytox PD-1+','T cytox Tim3-1+','T cytox PD1+Tim3+','Others GzmB/Perf+'};
% 
% count = 0;
% figflag = 1;
% 
% for i = 1:maxmouse
%     i
%     rois_ind = find(options.MouseNum==265+i);
%     index = ismember(MorpResults.Indexes,rois_ind) & ROIResults.SpleenIndex == 0;
%         
%     %number of tumors
%     tumor_id = unique(ROIResults.TumorIndex(index));
%     tumor_num = length(tumor_id);
%     for j = 1:tumor_num
%         temp = [];
%         totcells{i}(j) = sum(ROIResults.TumorIndex(index)==tumor_id(j));
%         epicells{i}(j) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,1)==2);
%         prolifer{i}(j) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,2)==21)/epicells{i}(j)*100;
%         
%         norm = totcells{i}(j);
%         immune{i,1}(j) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,3)==111)/norm*100; %t_cells
%         immune{i,2}(j) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,3)==112)/norm*100; % b_cells
%         immune{i,3}(j) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,3)==123)/norm*100; % nk_cells
%         immune{i,4}(j) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,3)==121)/norm*100; % macrophage
%         immune{i,5}(j)  = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,3)==122)/norm*100; % dendritic
%         immune{i,6}(j)  = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,2)==13)/norm*100; % neutrophil
%         immune_names = {'T','B','NK','Macro','Dendr','Neutro'};
%         
%         tcells{i,1}(j) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,4)==1111)/norm*100; %t_cells
%         tcells{i,2}(j) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,4)==1112)/norm*100; %t_cells
%         tcells{i,3}(j) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,4)==1113)/norm*100; %t_cells
%         tcells_names = {'T reg','T helper','T cytox'};
%         
%         temp(1) = i; % mouse number
%         temp(2) = totcells{i}(j); % total cell number
%         temp(3) = epicells{i}(j); % epith cell number
%         temp(4) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,3)==111); % T cell number
%         temp(5) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,3)==112); % B cell number
%         temp(6) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,3)==123); % nk_cells
%         temp(7) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,3)==121); % macrophage
%         temp(8) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,3)==122); % dendritic
%         temp(9) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,2)==13); % neutrophil        
%         temp(10) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,4)==1111); % treg mouse number
%         temp(11) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,4)==1112); % Thelper mouse number
%         temp(12) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,4)==1113); % Tcytotox mouse number
%         temp(13) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,4)==1113 & NormResults.MedianNucNorm(index,find(strcmp(label_mat{13}(9:end-1),options.Markers))) > 0); % Tcytotox mouse number
%         temp(14) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,4)==1113 & NormResults.MedianNucNorm(index,find(strcmp(label_mat{14}(9:end-1),options.Markers))) > 0); % Tcytotox mouse number
%         temp(15) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,4)==1113 & NormResults.MedianNucNorm(index,find(strcmp(label_mat{15}(9:end-1),options.Markers))) > 0); % Tcytotox mouse number
%         temp(16) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,4)==1113 & NormResults.MedianCytNorm(index,find(strcmp(label_mat{16}(9:end-1),options.Markers))) > 0); % Tcytotox mouse number
%         temp(17) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,4)==1113 & NormResults.MedianCytNorm(index,find(strcmp(label_mat{17}(9:end-1),options.Markers))) > 0); % Tcytotox mouse number
%         temp(18) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,4)==1113 & NormResults.MedianCytNorm(index,find(strcmp(label_mat{16}(9:end-1),options.Markers))) > 0 & NormResults.MedianCytNorm(index,find(strcmp(label_mat{17}(9:end-1),options.Markers))) > 0); % Tcytotox mouse number
% %         temp(19) = sum(ROIResults.TumorIndex(index)==tumor_id(j) & CellType.Matrix(index,2)==14); % Others GzmB/Perf+
%         
%         if tumor_id(j) ~= 0
%             count = count + 1;
%             mat(count,:) = temp;
%         else
%             nontumor(i,:) = temp;
%         end
%         
%     end
%     
%     prolifer_nontumor(i) = prolifer{i}(1);
%     prolifer_mean(i) = nanmedian(prolifer{i}(2:end));
%     prolifer_pr25(i) = prctile(prolifer{i}(2:end),25);
%     prolifer_pr75(i) = prctile(prolifer{i}(2:end),75);
%     
%     if figflag
%         figure(110)
%         bar(i,nanmedian(prolifer{i}(2:end)),colorlist{ceil(i/size(forplot.group_mat,2))})
%         hold on
%         plot([i i],[prctile(prolifer{i}(2:end),25) prctile(prolifer{i}(2:end),75)],colorlist{ceil(i/size(forplot.group_mat,2))},'LineWidth',3)
%         scatter(i,prolifer{i}(1),25*dotsizemult,'k','filled')
%         title('Proliferating cells per tumor')
% 
%         if mod(i,5)==1
%             text(i-0.5,-1,forplot.grouplist{ceil(i/size(forplot.group_mat,2))},'Color',colorlist{ceil(i/size(forplot.group_mat,2))})
%         end
%     end
%     for type = 1:6
%         label_all{type} = immune_names{type};
%         immune_nontumor(i,type) = immune{i,type}(1);
%         immune_mean(i,type) = nanmedian(immune{i,type}(2:end));
%         immune_pr25(i,type) = prctile(immune{i,type}(2:end),25);
%         immune_pr75(i,type) = prctile(immune{i,type}(2:end),75);
%         
%         if figflag
%             figure(120)
%             subplot(3,3,type)
%             xlim([0.5 maxmouse+0.5])
%             bar(i,nanmedian(immune{i,type}(2:end)),colorlist{ceil(i/size(forplot.group_mat,2))})
%             hold on
%             plot([i i],[prctile(immune{i,type}(2:end),25) prctile(immune{i,type}(2:end),75)],colorlist{ceil(i/size(forplot.group_mat,2))},'LineWidth',3)
%             scatter(i,immune{i,type}(1),25*dotsizemult,'k','filled')
%             if i == maxmouse
%                 title(immune_names{type})
%             end
%             if mod(i,size(forplot.group_mat,2))==3
%                 text(i/20,-0.1,forplot.grouplist{ceil(i/size(forplot.group_mat,2))},'Color',colorlist{ceil(i/size(forplot.group_mat,2))}, ...
%                      'Units','normalized','HorizontalAlignment','center')
%             end
%             if type == 4
%                 ylabel('% cells of total cell in tumor area (%) - median +/- 25th percentile')
%             end
%         end
%         
%     end
%     for type = 1:3
%         label_all{6+type} = tcells_names{type};
%         immune_nontumor(i,6+type) = tcells{i,type}(1);
%         immune_mean(i,6+type) = nanmedian(tcells{i,type}(2:end));
%         immune_pr25(i,6+type) = prctile(tcells{i,type}(2:end),25);
%         immune_pr75(i,6+type) = prctile(tcells{i,type}(2:end),75);
%         
%         if figflag
%             figure(120)
%             subplot(3,3,type+6)
%             xlim([0.5 maxmouse+0.5])
%             bar(i,nanmedian(tcells{i,type}(2:end)),colorlist{ceil(i/size(forplot.group_mat,2))})
%             hold on
%             plot([i i],[prctile(tcells{i,type}(2:end),25) prctile(tcells{i,type}(2:end),75)],colorlist{ceil(i/size(forplot.group_mat,2))},'LineWidth',3)
%             scatter(i,tcells{i,type}(1),25*dotsizemult,'k','filled')
%             if mod(i,size(forplot.group_mat,2))==3
%                 text(i/20,-0.1,forplot.grouplist{ceil(i/size(forplot.group_mat,2))},'Color',colorlist{ceil(i/size(forplot.group_mat,2))}, ...
%                      'Units','normalized','HorizontalAlignment','center')
%             end
%             if i == maxmouse
%                 title(tcells_names{type})
%             end
%         end
%     end
% end
% figure(120)
% set(gcf,'Units','normalized','Position',[0 0.05 1 0.85])
% %%
% % compare statistically the cell population numbers between Cre dgTom/dgCxcl10
% % between Cre and LucOS
% 
% size_th = 10;
% rep = [];
% for type = 4:18
%     rep{type,1} = label_mat{type};
%     for i = 1:4
%         
%         ind = ceil(mat(:,1)/5) == i & mat(:,2) > size_th;
%         k{i} = mat(ind,type)./mat(ind,3);
%         
%         rep{type,i+1} = median(k{i}(k{i}>0))*100;
%         
%         ind = ceil(nontumor(:,1)/5) == i;
%         rep{type,7+i} = median(nontumor(ind,type)./nontumor(ind,3))*100;
%         
%     end
%     
%     [a,b]=kstest2(k{1}(k{1}>0),k{2}(k{2}>0));
%     rep{type,6} = b;
%     [a,b]=kstest2(k{1}(k{1}>0),k{3}(k{3}>0));
%     rep{type,7} = b;
%     
% end
% 
% 
% %%
% colors_list = {'c','b','m','r'};
% figure(300)
% for i = [1 2]
%     s = (i-1)*5+1;
%     e = s+4;
%     scatter(prolifer_nontumor(s:e),prolifer_mean(s:e),45,colors_list{i},'filled')
%     hold on
% end
% plot([0 max([prolifer_nontumor prolifer_mean])*1.05],[0 max([prolifer_nontumor prolifer_mean])*1.05],'--k')
% axis([0 max([prolifer_nontumor prolifer_mean])*1.05 0 max([prolifer_nontumor prolifer_mean])*1.05])
% title('% proliferating cells')
% legend('Cre','dgCxcl10','LucOS')
% ylabel('Cell Frequency INSIDE tumor areas')
% xlabel('Cell Frequency OUTSIDE tumor areas')
% 
% for j = 1:9
%     figure(310)
%     subplot(3,3,j)
%     
%     for i = [1 2]
%         s = (i-1)*5+1;
%         e = s+4;
%         scatter(immune_nontumor(s:e,j),immune_mean(s:e,j),45,colors_list{i},'filled')
%         hold on
%     end
%     plot([0 max([immune_nontumor(:,j); immune_mean(:,j)])*1.05],[0 max([immune_nontumor(:,j); immune_mean(:,j)])*1.05],'--k')
%     axis([0 max([immune_nontumor(:,j); immune_mean(:,j)])*1.05 0 max([immune_nontumor(:,j); immune_mean(:,j)])*1.05])
%     title(label_all{j})
%     
%     if j == 4
%         ylabel('Cell Frequency INSIDE tumor areas')
%     elseif j == 8
%         xlabel('Cell Frequency OUTSIDE tumor areas')
%     elseif j == 9
%         legend('Cre','LucOS')
%     end
%     
%     figure(320)
%     subplot(3,3,j)
%     for i = [1 2]
%         s = (i-1)*5+1;
%         e = s+4;
%         scatter(immune_nontumor(s:e,j),immune_mean(s:e,j),45,colors_list{i},'filled')
%         hold on
%     end
%     plot([0 max([immune_nontumor(:,j); immune_mean(:,j)])*1.05],[0 max([immune_nontumor(:,j); immune_mean(:,j)])*1.05],'--k')
%     axis([0 max([immune_nontumor(:,j); immune_mean(:,j)])*1.05 0 max([immune_nontumor(:,j); immune_mean(:,j)])*1.05])
%     title(label_all{j})
%     
%     if j == 4
%         ylabel('Cell Frequency INSIDE tumor areas')
%     elseif j == 8
%         xlabel('Cell Frequency OUTSIDE tumor areas')
%     elseif j == 9
%         legend('dgTom','dgCxcl10')
%     end
% end
% figure(310)
% set(gcf,'Units','normalized','Position',[0.2 0.05 0.4 0.8])
% figure(320)
% set(gcf,'Units','normalized','Position',[0.2 0.05 0.4 0.8])


% %% basic stats regarding networks
% clear NetsMouse
% plotvect = [ones(5,1)-(6-(1:5)'-3)*0.1; 2*ones(5,1)-(6-(1:5)'-3)*0.1];
% net_thresh = 20;
% for i = 1:10
%     rois_ind = find(options.MouseNum==forplot.group_mat(ceil(i/5),i-(ceil(i/5)-1)*5));
%     
%     NetsMouse.NumCells{i} = [];   
%     NetsMouse.Bfrac{i} = [];
%     NetsMouse.cytoTfrac{i} = [];
%     NetsMouse.regTfrac{i} = [];
%     NetsMouse.helpTfrac{i} = [];
%     NetsMouse.otherfrac{i} = [];
%     
%     NetsMouse.NetsNum(i) = 0;
%     NetsMouse.CellInNets(i,1:3) = 0;
%     temp_totcells = 0;
%     for r = 1:length(rois_ind)
%         try 
%             ind_net = LymphoNets.Summary{rois_ind(r)}.NumCells > net_thresh;
%         catch
%             continue
%         end
%         NetsMouse.NumCells{i} = [NetsMouse.NumCells{i}; LymphoNets.Summary{rois_ind(r)}.NumCells(ind_net)'];
%         NetsMouse.NetsNum(i) = NetsMouse.NetsNum(i) + sum(ind_net);
%         
%         % extract total number of lymphocytes
%         lymphocytes = sum(MorpResults.Indexes == rois_ind(r) & ROIResults.SpleenIndex == 0 & CellType.Matrix(:,2) == 11);
%         totcells = sum(MorpResults.Indexes == rois_ind(r) & ROIResults.SpleenIndex == 0);
%         NetsMouse.CellInNets(i,2) = NetsMouse.CellInNets(i,2) + lymphocytes;
%         temp_totcells = temp_totcells + totcells;
%         
%         % B vs T cells
%         NetsMouse.Bfrac{i} =     [NetsMouse.Bfrac{i};     LymphoNets.Summary{rois_ind(r)}.Bcells(ind_net)'./LymphoNets.Summary{rois_ind(r)}.NumCells(ind_net)'];
%         NetsMouse.cytoTfrac{i} = [NetsMouse.cytoTfrac{i}; LymphoNets.Summary{rois_ind(r)}.Tcytox(ind_net)'./LymphoNets.Summary{rois_ind(r)}.NumCells(ind_net)'];
%         NetsMouse.regTfrac{i} =  [NetsMouse.regTfrac{i}; LymphoNets.Summary{rois_ind(r)}.Tregs(ind_net)' ./LymphoNets.Summary{rois_ind(r)}.NumCells(ind_net)'];
%         NetsMouse.helpTfrac{i} = [NetsMouse.helpTfrac{i}; LymphoNets.Summary{rois_ind(r)}.Thelps(ind_net)'./LymphoNets.Summary{rois_ind(r)}.NumCells(ind_net)'];
%         NetsMouse.otherfrac{i} = [NetsMouse.otherfrac{i}; 1 - NetsMouse.Bfrac{i} - NetsMouse.cytoTfrac{i} - NetsMouse.regTfrac{i} - NetsMouse.helpTfrac{i} ];
%     end
%     NetsMouse.CellInNets(i,1) = sum(NetsMouse.NumCells{i});
%     NetsMouse.CellInNets(i,3) = NetsMouse.CellInNets(i,1)./NetsMouse.CellInNets(i,2);
%     NetsMouse.CellInNets(i,4) = NetsMouse.CellInNets(i,2)./temp_totcells;
%     
%     NetsMouse.NetSize(i,1) = mean(NetsMouse.NumCells{i});
%     NetsMouse.NetSize(i,2) = prctile(NetsMouse.NumCells{i},25);
%     NetsMouse.NetSize(i,3) = prctile(NetsMouse.NumCells{i},75);
%     
%     NetsMouse.BfracSummary(i,1) = mean(NetsMouse.Bfrac{i});
%     NetsMouse.BfracSummary(i,2) = prctile(NetsMouse.Bfrac{i},25);
%     NetsMouse.BfracSummary(i,3) = prctile(NetsMouse.Bfrac{i},75);
%     
%     NetsMouse.cytoTfracSummary(i,1) = mean(NetsMouse.cytoTfrac{i});
%     NetsMouse.regTfracSummary(i,1) = mean(NetsMouse.regTfrac{i});
%     NetsMouse.helpTfracSummary(i,1) = mean(NetsMouse.helpTfrac{i});
%     NetsMouse.otherfracSummary(i,1) = mean(NetsMouse.otherfrac{i});
%     
% %     figure(3000)
% %     [n,h] = hist(NetsMouse.NumCells{i},((1:21)-1)*5);
% %     subplot(2,2,ceil(i/5))
% %     plot(h,n/sum(n))
% %     hold on
%     
% end
% 
% for m = 1:2
%     miceintreat = ((m-1)*6+1):(m*6);
%     mean_forbar_numnets(m) = nanmean(NetsMouse.NetsNum(miceintreat));
%     mean_forbar_lymphnum(m) = nanmean(NetsMouse.CellInNets(miceintreat,4));
%     mean_forbar_lymphfrac(m) = nanmean(NetsMouse.CellInNets(miceintreat,3));
%     mean_forbar_netsize(m) = nanmean(NetsMouse.NetSize(miceintreat,1));
%     mean_forbar_bfrac(m) = nanmean(NetsMouse.BfracSummary(miceintreat,1));
%     mean_forbar_cytotfrac(m) = nanmean(NetsMouse.cytoTfracSummary(miceintreat,1));
%     mean_forbar_regtfrac(m) = nanmean(NetsMouse.regTfracSummary(miceintreat,1));
%     mean_forbar_helptfrac(m) = nanmean(NetsMouse.helpTfracSummary(miceintreat,1));
% end
% 
%        
% figure
% subplot(2,4,1)
% bar([1:length(mean_forbar_lymphnum)],mean_forbar_numnets,0.5,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
% hold on
% scatter(plotvect,NetsMouse.NetsNum,40,round(plotvect),'filled')
% caxis([0 4.5])
% xlim([0.5 4.5])
% ylim([0 max(NetsMouse.NetsNum)*1.1])
% colormap('jet')
% set(gca,'XTick',1:4,'XTickLabel',forplot.grouplist,'XTickLabelRotation',45)
% title('Number of lymphonets')
% ylabel('# networks in 2 lobes')
% 
% subplot(2,4,2)
% bar([1:length(mean_forbar_lymphnum)],mean_forbar_lymphnum,0.5,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
% hold on
% scatter(plotvect,NetsMouse.CellInNets(:,4),40,round(plotvect),'filled')
% caxis([0 4.5])
% xlim([0.5 4.5])
% ylim([0 max(NetsMouse.CellInNets(:,4))*1.1])
% colormap('jet')
% set(gca,'XTick',1:4,'XTickLabel',forplot.grouplist,'XTickLabelRotation',45)
% title('Total lymph pop')
% ylabel('# lymphocytes/ # tot cells')
% 
% subplot(2,4,3)
% bar([1:length(mean_forbar_lymphnum)],mean_forbar_lymphfrac,0.5,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
% hold on
% scatter(plotvect,NetsMouse.CellInNets(:,3),40,round(plotvect),'filled')
% caxis([0 4.5])
% xlim([0.5 4.5])
% ylim([0 max(NetsMouse.CellInNets(:,3))*1.1])
% colormap('jet')
% set(gca,'XTick',1:4,'XTickLabel',forplot.grouplist,'XTickLabelRotation',45)
% title('Frac lymph in Nets')
% ylabel('# in nets / tot')
% 
% subplot(2,4,4)
% bar([1:length(mean_forbar_lymphnum)],mean_forbar_netsize,0.5,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
% hold on
% scatter(plotvect,NetsMouse.NetSize(:,1),40,round(plotvect),'filled')
% caxis([0 4.5])
% xlim([0.5 4.5])
% ylim([0 max(NetsMouse.NetSize(:,1))*1.1])
% colormap('jet')
% set(gca,'XTick',1:4,'XTickLabel',forplot.grouplist,'XTickLabelRotation',45)
% title('Mean Net Size')
% ylabel('# cells')
% 
% subplot(2,4,5)
% bar([1:length(mean_forbar_lymphnum)],mean_forbar_bfrac,0.5,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
% hold on
% scatter(plotvect,NetsMouse.BfracSummary(:,1),40,round(plotvect),'filled')
% caxis([0 4.5])
% xlim([0.5 4.5])
% ylim([0 max(NetsMouse.BfracSummary(:,1))*1.1])
% colormap('jet')
% set(gca,'XTick',1:4,'XTickLabel',forplot.grouplist,'XTickLabelRotation',45)
% title('Mean B cell frac')
% 
% subplot(2,4,6)
% bar([1:length(mean_forbar_lymphnum)],mean_forbar_cytotfrac,0.5,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
% hold on
% scatter(plotvect,NetsMouse.cytoTfracSummary(:,1),40,round(plotvect),'filled')
% caxis([0 4.5])
% xlim([0.5 4.5])
% ylim([0 max(NetsMouse.cytoTfracSummary(:,1))*1.1])
% colormap('jet')
% set(gca,'XTick',1:4,'XTickLabel',forplot.grouplist,'XTickLabelRotation',45)
% title('Mean cytoT cell frac')
% 
% subplot(2,4,7)
% bar([1:length(mean_forbar_lymphnum)],mean_forbar_regtfrac,0.5,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
% hold on
% scatter(plotvect,NetsMouse.regTfracSummary(:,1),40,round(plotvect),'filled')
% caxis([0 4.5])
% xlim([0.5 4.5])
% ylim([0 max(NetsMouse.regTfracSummary(:,1))*1.1])
% colormap('jet')
% set(gca,'XTick',1:4,'XTickLabel',forplot.grouplist,'XTickLabelRotation',45)
% title('Mean regT cell frac')
% 
% subplot(2,4,8)
% bar([1:length(mean_forbar_lymphnum)],mean_forbar_helptfrac,0.5,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
% hold on
% scatter(plotvect,NetsMouse.helpTfracSummary(:,1),40,round(plotvect),'filled')
% caxis([0 4.5])
% xlim([0.5 4.5])
% ylim([0 max(NetsMouse.helpTfracSummary(:,1))*1.1])
% colormap('jet')
% set(gca,'XTick',1:4,'XTickLabel',forplot.grouplist,'XTickLabelRotation',45)
% title('Mean helpT cell frac')
% 
%% calculate the distance of cells from the border of the closest tumor and compare it to the distance from the closest blood vessel

% strategy:
% - load mask from tumors and blood vessels
% 1) measure the distance from blood vessel
% 2) measure the distance from a tumor
% 3) measure the distance from the border of a tumor

DistResults.Tumor = zeros(size(MorpResults.Indexes,1),3); 
% col (1) tumor #
% col (2) distance from tumor - 0 if inside
% col (3) distance from tumor boundary
% DistResults.Blood = zeros(size(MorpResults.Indexes,1),1); 
% col (1) distance from blood vessel - 0 if inside
% col (2) distance from tumor boundary

for tis = 1:max(MorpResults.Indexes)
    tis
    % load masks for tumor and blood
    basename = [filename.analfolder filename.roifolder filename.folders{tis} '_'];
    TumorMask = imread([basename options.PathLabels{1} '.tif']);
%     BloodMask = imread([basename options.tifnames.Path{2}]);
    TumorBoundMask = bwperim(TumorMask);
%     BloodBoundMask = bwperim(BloodMask);
    
    % calcuate the distance
    [TumorDist, idx] = bwdist(TumorMask);
%     BloodDist = bwdist(BloodMask);
    TumorBoundDist = bwdist(TumorBoundMask);
%     BloodBoundDist = bwdist(BloodBoundMask);
    
    
    % get cell index
    vect_cell_ind = MorpResults.Indexes==tis;
    X_scale = max([MorpResults.X(vect_cell_ind)/(2^options.ROI.pyramidlevel) ones(sum(vect_cell_ind),1)],[],2);
    Y_scale = max([MorpResults.Y(vect_cell_ind)/(2^options.ROI.pyramidlevel) ones(sum(vect_cell_ind),1)],[],2);
    cell_index = sub2ind(size(TumorMask),Y_scale,X_scale);
    
%     figure
%     subplot(1,2,1)
%     imshow(TumorMask,[])
%     subplot(1,2,2)
%     scatter(X_scale,Y_scale,4,'filled')
%     set(gca,'YDir','reverse')
%     axis([0 size(TumorMask,2) 0 size(TumorMask,1)])
    
    % calculate distance matrices
    DistResults.Tumor(vect_cell_ind,1) = TumorMask(idx(cell_index));
    DistResults.Tumor(vect_cell_ind,2) = TumorDist(cell_index);
    DistResults.Tumor(vect_cell_ind,3) = TumorBoundDist(cell_index);
    
%     DistResults.Blood(vect_cell_ind,1) = BloodDist(cell_index);
%     DistResults.Blood(vect_cell_ind,2) = BloodBoundDist(cell_index);
    
    % now we do it for the centers of the networks of lymphocytes
    % initialize matrices for networks of lymphocytes
    LymphoNets.Summary{tis}.DistTumor =  double(zeros(size(LymphoNets.Summary{tis}.X_mean,2),3)); 
%     LymphoNets.Summary{tis}.DistBlood =  double(zeros(size(LymphoNets.Summary{tis}.X_mean,2),2)); 
    
    % get cell index
    X_ln_scale = round(LymphoNets.Summary{tis}.X_mean'/(2^options.ROI.pyramidlevel),0);
    Y_ln_scale = round(LymphoNets.Summary{tis}.Y_mean'/(2^options.ROI.pyramidlevel),0);
    ln_index = sub2ind(size(TumorMask),Y_ln_scale,X_ln_scale);
    
    LymphoNets.Summary{tis}.DistTumor(:,1) = TumorMask(idx(ln_index));
    LymphoNets.Summary{tis}.DistTumor(:,2) = TumorDist(ln_index);
    LymphoNets.Summary{tis}.DistTumor(:,3) = TumorBoundDist(ln_index);
%     LymphoNets.Summary{tis}.DistBlood(:,1) = BloodDist(ln_index);
%     LymphoNets.Summary{tis}.DistBlood(:,2) = BloodBoundDist(ln_index);
    
    LymphoNets.Summary{tis}.DistTumor(LymphoNets.Summary{tis}.DistTumor(:,2)==0,3) = - LymphoNets.Summary{tis}.DistTumor(LymphoNets.Summary{tis}.DistTumor(:,2)==0,3); 
%     LymphoNets.Summary{tis}.DistBlood(LymphoNets.Summary{tis}.DistBlood(:,1)==0,2) = - LymphoNets.Summary{tis}.DistBlood(LymphoNets.Summary{tis}.DistBlood(:,1)==0,2);
end
% make the distance from boundary negative if the cell is inside the tumor
DistResults.Tumor(DistResults.Tumor(:,2)==0,3) = -DistResults.Tumor(DistResults.Tumor(:,2)==0,3);
% DistResults.Blood(DistResults.Blood(:,1)==0,2) = -DistResults.Blood(DistResults.Blood(:,1)==0,2);

%%
% %% plot all the subtype distances
% umperpix = 0.325;
% for i = 1:20
%     i
%     rois_ind = find(options.MouseNum==265+i);
%     tissue = ismember(MorpResults.Indexes,rois_ind) & ROIResults.SpleenIndex == 0;
%     for type = 1:length(CellType.codes)
%         % take all the tcells
%         layer_temp = ceil(log10(double(CellType.codes(type))+0.1));
%         subtype = CellType.Matrix(:,layer_temp) == CellType.codes(type);
%         subtype_dist = DistResults.Tumor(tissue & subtype,3)*umperpix;
%         
%         figure(10000+type)
%         subplot(2,2,ceil(i/5))
%         [n,h] = hist(subtype_dist,linspace(-50,120,18));
%         plot(h,n/sum(n))
%         hold on
%         title([CellType.names{type} ' ' forplot.grouplist{ceil(i/5)}])
%         xlabel('distance from tumor border (um)')
%         axis([-55 125 0 0.65])
%     end
% end
% 
% %% Phenotyping of T cell w respect to position to tumor
% 
% marker_list = {'GranzB','PD-1','Tim3-2','Ki67'};
% thresh =      [ 0      0    0      0];%[ 300      200    200      200];
% celltype_list = {'T cytotox'};%,'T reg','T helper'};
% treat_list = {'Luc','LucOS'};
% 
% for i1 = 1:length(marker_list)
%     for i2 = 1:length(celltype_list)
%         
%         % find the right cell type
%         type_pos = find(strcmp(CellType.names,celltype_list{i2}));
%         type_code = CellType.codes(type_pos);
%         type_layer = CellType.layer(type_pos);
%         cells = CellType.Matrix(:,type_layer) == type_code;
% 
%         % find the right marker
%         mark = find(strcmp(options.Markers,marker_list{i1}));
%         figure(4000)
%         
%         for i3 = [1 3]
%             mouse_group = find(options.MouseGroup==i3);
% 
%             ind = ismember(MorpResults.Indexes,mouse_group) & cells & ROIResults.SpleenIndex == 0;
%             v_x = NormResults.MedianNucNorm(ind,mark);
%             v_y = DistResults.Tumor(ind,3);
%             v_y = log10(double(abs(v_y)+1)).*sign(v_y);
%             bin_dist = linspace(-2,3,21);
%             
%             subplot(length(marker_list),2*length(celltype_list),(i1-1)*length(celltype_list)*length(treat_list)+(i2-1)*length(treat_list)+ceil(i3/2))
%             [n,h]=hist(v_y,bin_dist);
%             plot(h,n/sum(n))
%             hold on
%             [n,h]=hist(v_y(v_x>thresh(i1)),bin_dist);
%             plot(h,n/sum(n))
%             plot([0 0],[0 0.2],'--k',[1.2 1.2],[0 0.2],'--k')
%             text(-1.6,0.22,'Inside')
%             text(0.07,0.22,'Prox')
%             text(2.5,0.22,'Distal')
%             text(3.2,0.18,['n=' num2str(length(v_x))])
%             text(3.2,0.12,[ num2str(round(sum(v_x>thresh(i1))/length(v_x)*100,2)) '% +ve'])
%             title([celltype_list{i2} ' ' marker_list{i1} ' ' treat_list{ceil(i3/2)}])
%             set(gcf,'Units','Normalized','OuterPosition',[0 0.04 1 0.96])
%             axis([-2 6 0 0.24])
%             
%         end
% 
%     end
% end
% 
% %% plot of lymphocyte networks with respect to tumors and blood vessels
% colorbygroup = {'b','c','r','m','g'};
% PixInMicron = 2048/665; 
% net_thresh = 15;
% for j = 1:28
%     % first we need to plot the blood vessels and the tumors
%     tu = imread([filename.analfolder filename.roifolder filename.folders{j}(5:end) '_' options.tifnames.Path{1}]);
%     bl = imread([filename.analfolder filename.roifolder filename.folders{j}(5:end) '_' options.tifnames.Path{2}]);
%    
%     % use all the non-spleen points to make a image
%     index = MorpResults.Indexes == j & ROIResults.SpleenIndex == 0;
%     Y_scale = round(MorpResults.Y(index)/2^options.ROI.pyramidlevel);
%     X_scale = round(MorpResults.X(index)/2^options.ROI.pyramidlevel);
%     whole = uint16(zeros(size(tu)));
%     graph_index = sub2ind(size(whole),Y_scale,X_scale);
%     whole(graph_index)=1;
%     whole = imdilate(whole,strel('disk',3));
%     whole = imfill(whole,'holes');
%     whole = imerode(whole,strel('disk',3));
%     
%     % combine the matrices so that the blood vessels show up on top of
%     % tumors
%     comp = cat(3,whole,uint16(tu>0)*2,uint16(bl*3));
%     comp = max(comp,[],3);
%     R = imref2d(size(tu),[0 size(tu,2)*2^(options.ROI.pyramidlevel)],[0 size(tu,1)*2^(options.ROI.pyramidlevel)]);
%     
% %     figure(5000+j)
% %     imshow(comp,R,[0 4]); colormap(flipud(gray));
% %     hold on
% %     scatter(MorpResults.X(index & LymphoNets.NetworkID > 0),MorpResults.Y(index  & LymphoNets.NetworkID > 0),20,[255 128 0]/255,'filled')
% %     viscircles([LymphoNets.Summary{j}.X_mean' LymphoNets.Summary{j}.Y_mean'],LymphoNets.Summary{j}.NumCells','Color','b');
% %     title(filename.tissues{j})
%     
%     ind_net = LymphoNets.Summary{j}.NumCells > net_thresh;
%     
%     logdist_tumor = sign(LymphoNets.Summary{j}.DistTumor(ind_net,3)).*log2(abs(LymphoNets.Summary{j}.DistTumor(ind_net,3))+1);
%     logdist_blood = sign(LymphoNets.Summary{j}.DistBlood(ind_net,2)).*log2(abs(LymphoNets.Summary{j}.DistBlood(ind_net,2))+1);
%     n = 3;
%     lognumcell = (double(LymphoNets.Summary{j}.NumCells(ind_net)))/n;
%     
%     figure(6000+options.MouseGroup(j))
%     scatter(logdist_tumor,logdist_blood,lognumcell,colorbygroup{options.MouseGroup(j)},'filled','MarkerEdgeColor',[0.8 0.8 0.8])
%     hold on
%     plot([0 0],[-4 10],'-k')
%     plot([-8 10],[0 0],'-k')
%     plot([log2(PixInMicron*10) log2(PixInMicron*10)],[-4 10],'--r')
%     plot([-8 10],[log2(PixInMicron*10) log2(PixInMicron*10)],'--r')
%     if mod(options.MouseNum(j),5) == 0
%         title(forplot.grouplist{options.MouseGroup(j)})
%         axis([-8 10 -4 11])
%         text(-6,10.3,'Inside Tumor')
%         text(1,10.3,'Within 10uM')
%         text(5.5,10.3,'Distal from Tumor')
%         xlabel('log2 distance from closest tumor')
%         ylabel('log2 distance from closest blood vessel')
%         
%         legx = zeros(6,1)-6;
% %         legy = [-0.9 -1.1 -1.4 -1.7 -2.1 -2.6 -3.3];
%         legy = -(1:6)*0.55;
%         legz = (2.^[-5 -4 -3 -2 -1 0])*250/n;
%         
%         scatter(legx,legy,legz,colorbygroup{options.MouseGroup(j)},'filled','MarkerEdgeColor',[0.8 0.8 0.8])
%         text(legx'+0.4,legy',num2str(round(legz'*n,0)))
%     end
%     
% end

%% INVESTIGATING THE SOLIDITY OF THE NETWORKS

bcell_num = [];
tcell_num = [];
count = 0;

for j = 1:length(options.MouseNum)    
    j
    prop = [];
    list_nets = unique(LymphoNets.NetworkID(MorpResults.Indexes ==j));
    
    all_cells = MorpResults.Indexes == j & ROIResults.SpleenIndex == 0;
    tum_cells = MorpResults.Indexes == j & CellType.Matrix(:,2) ~= 11 & ROIResults.TumorIndex > 0 ;
    
    coords = double([MorpResults.X(all_cells) MorpResults.Y(all_cells)]);    
    % perform delaunay graph for all the cells
    
%     if j == 15
%         figure(j)
%         scatter(MorpResults.X(all_cells),MorpResults.Y(all_cells),35,[0.8 0.8 0.8],'filled')
%         hold on
%         scatter(MorpResults.X(tum_cells),MorpResults.Y(tum_cells),35,[0.6 0.6 0.6],'filled')
%     end 
    
    for n = 2:length(list_nets)
        cell_ind = MorpResults.Indexes ==j & LymphoNets.NetworkID == list_nets(n);
   
        if sum(cell_ind) > 5
            count = count + 1;
            cum_solid(count,1) = options.MouseGroup(j);
            cum_solid(count,2) = sum(cell_ind);
%             cum_solid(count,3) = mean( double(LymphoNets.LymphoEdges(cell_ind))./ double(LymphoNets.AllEdges(cell_ind)));
            cum_solid(count,3) = (double(sum(LymphoNets.LymphoEdges(cell_ind)))/2) / (sum(double(LymphoNets.AllEdges(cell_ind))) - sum(double(LymphoNets.LymphoEdges(cell_ind)))/2 );
            
%             if j == 15
%                 figure(j)
%                 scatter(MorpResults.X(cell_ind),MorpResults.Y(cell_ind),35,cum_solid(count,3)+zeros(sum(cell_ind),1),'filled')
%                 caxis([0.4 1])
%                 colormap('cool')
%             end
            
            cum_ki67(count,1) = median(NormResults.MedianNucNorm(cell_ind,find(strcmp(options.Markers,'Ki67'))));
            cum_ki67(count,2) = prctile(NormResults.MedianNucNorm(cell_ind,find(strcmp(options.Markers,'Ki67'))),80);
            
            dist_tumor(count,1) = min(DistResults.Tumor(cell_ind,3));
%             dist_blood(count,1) = min(DistResults.Blood(cell_ind,2));
            
%             t_cell_ind = cell_ind & CellType.Matrix(:,3) == 111;
%             
%             cum_pd1(count,1) = median(NormResults.MedianNucNorm(t_cell_ind,find(strcmp(options.Markers,'PD-1'))));
%             cum_pd1(count,2) = prctile(NormResults.MedianNucNorm(t_cell_ind,find(strcmp(options.Markers,'PD-1'))),80);
            
%             cum_perf(count,1) = median(NormResults.MedianNucNorm(t_cell_ind,find(strcmp(options.Markers,'Perforin'))));
%             cum_perf(count,2) = prctile(NormResults.MedianNucNorm(t_cell_ind,find(strcmp(options.Markers,'Perforin'))),80);
            
%             cum_tim3(count,1) = median(NormResults.MedianNucNorm(t_cell_ind,find(strcmp(options.Markers,'Tim3-1'))));
%             cum_tim3(count,2) = prctile(NormResults.MedianNucNorm(t_cell_ind,find(strcmp(options.Markers,'Tim3-1'))),80);
            
            
            bcell_num(count) = sum(CellType.Matrix(cell_ind,3) == 112);
            tcell_num(count) = sum(CellType.Matrix(cell_ind,3) == 111);
            treg_num(count) = sum(CellType.Matrix(cell_ind,4) == 1111);
            thelp_num(count) = sum(CellType.Matrix(cell_ind,4) == 1112);
            tcyto_num(count) = sum(CellType.Matrix(cell_ind,4) == 1113);
            
            check =  sum(cell_ind) - sum(CellType.Matrix(cell_ind,2)==11);
            if check ~= 0
                disp('Something is not right')
                break
            end
        end
    end

end

%% 

jumps = [ 5 15 40 200 ];

for i = [1 2 3] 
    if i ~= 3
        ind = cum_solid(:,1) == i;
        title_for_plot = forplot.grouplist{i};
    else
        ind = cum_solid(:,1) < i;
        title_for_plot = 'All Cells';
    end
        
    for j = 1:length(jumps)-1
            ind_sub = ind & cum_solid(:,2) >= jumps(j) & cum_solid(:,2) < jumps(j+1) ;
            if sum(ind_sub) > 15 
                figure(20)
                subplot(2,2,i)
                scatter(log2(double(cum_solid(ind_sub,2))),cum_solid(ind_sub,3),10,[0.8 0.8 0.8],'filled'); 
                hold on
                [f,x] = ksdensity([log2(double(cum_solid(ind_sub,2))) cum_solid(ind_sub,3)]);%,'PlotFcn','contour');
                f = f/sum(f);
                [M,c] = contour(reshape(x(:,1),[30 30]),reshape(x(:,2),[30 30]),reshape(f,[30 30]),40);
                colormap('jet')
                c.LineWidth = 2;
                ylim([0.2 1])
                xlim([2 8])
                xlabel('log2 number of cells')
                ylabel('Network Solidity')
                title(title_for_plot )
            end
%             
%             % separate cells by PD-1 top quintile ONLY OF T CELLS!
%             ind_sub = ind & cum_pd1(:,2) < 0  & cum_solid(:,2) >= jumps(j) & cum_solid(:,2) < jumps(j+1) ; % cum_ki67 cum_pd1(:,2) < 0   dist_tumor < 0 
%             figure(200)
%             subplot(4,2,1+(i-1)*2)
%             scatter(log2(double(cum_solid(ind_sub,2))),cum_solid(ind_sub,3),10,[0.8 0.8 0.8],'filled'); 
%             hold on
%             if sum(ind_sub) > 15 
%                 [f,x] = ksdensity([log2(double(cum_solid(ind_sub,2))) cum_solid(ind_sub,3)]);%,'PlotFcn','contour');
%                 [M,c] = contour(reshape(x(:,1),[30 30]),reshape(x(:,2),[30 30]),reshape(f,[30 30]),20);
%                 colormap('jet')
%                 c.LineWidth = 2;
%             end
%             ylim([0.2 1])
%             xlim([2 8])
%             title([title_for_plot ' T-cell PD-1 top quintile < 0'])
%             xlabel('log2 cell num')
% 
%             ind_sub = ind  & cum_pd1(:,2) >= 0 & cum_solid(:,2) >= jumps(j) & cum_solid(:,2) < jumps(j+1) ; %  cum_ki67 cum_pd1(:,2) >= 0 dist_tumor >= 0
%             figure(200)
%             subplot(4,2,2+(i-1)*2)
%             scatter(log2(double(cum_solid(ind_sub,2))),cum_solid(ind_sub,3),10,[0.8 0.8 0.8],'filled'); 
%             hold on
%             if sum(ind_sub) > 15 
%                 [f,x] = ksdensity([log2(double(cum_solid(ind_sub,2))) cum_solid(ind_sub,3)]);%,'PlotFcn','contour');
%                 [M,c] = contour(reshape(x(:,1),[30 30]),reshape(x(:,2),[30 30]),reshape(f,[30 30]),20);
%                 colormap('jet')
%                 c.LineWidth = 2;
%             end          
%             ylim([0.2 1])
%             xlim([2 8])
%             title([title_for_plot ' T-cell PD-1 top quintile > 0'])
%             xlabel('log2 cell num')
%             
%             % separate cells by Ki67 top quintile
%             ind_sub = ind & cum_ki67(:,2) < 0  & cum_solid(:,2) >= jumps(j) & cum_solid(:,2) < jumps(j+1) ; % cum_ki67 cum_perf cum_tim3
%             figure(300)
%             subplot(4,2,1+(i-1)*2)
%             scatter(log2(double(cum_solid(ind_sub,2))),cum_solid(ind_sub,3),10,[0.8 0.8 0.8],'filled'); 
%             hold on
%             if sum(ind_sub) > 15 
%                 [f,x] = ksdensity([log2(double(cum_solid(ind_sub,2))) cum_solid(ind_sub,3)]);%,'PlotFcn','contour');
%                 [M,c] = contour(reshape(x(:,1),[30 30]),reshape(x(:,2),[30 30]),reshape(f,[30 30]),20);
%                 colormap('jet')
%                 c.LineWidth = 2;
%             end
%             ylim([0.2 1])
%             xlim([2 8])
%             title([title_for_plot ' Ki67 top quintile < 0'])
%             xlabel('log2 cell num')
% 
%             ind_sub = ind  & cum_ki67(:,2) >= 0 & cum_solid(:,2) >= jumps(j) & cum_solid(:,2) < jumps(j+1) ; %  cum_ki67 cum_perf cum_tim3
%             figure(300)
%             subplot(4,2,2+(i-1)*2)
%             scatter(log2(double(cum_solid(ind_sub,2))),cum_solid(ind_sub,3),10,[0.8 0.8 0.8],'filled'); 
%             hold on
%             if sum(ind_sub) > 15 
%                 [f,x] = ksdensity([log2(double(cum_solid(ind_sub,2))) cum_solid(ind_sub,3)]);%,'PlotFcn','contour');
%                 [M,c] = contour(reshape(x(:,1),[30 30]),reshape(x(:,2),[30 30]),reshape(f,[30 30]),20);
%                 colormap('jet')
%                 c.LineWidth = 2;
%             end
%             ylim([0.2 1])
%             xlim([2 8])
%             title([title_for_plot ' Ki67 top quintile > 0'])
%             xlabel('log2 cell num')
%             
            % separate cells by location wrt tumor
            ind_sub = ind & dist_tumor < 0  & cum_solid(:,2) >= jumps(j) & cum_solid(:,2) < jumps(j+1) ; % cum_ki67 cum_pd1(:,2) < 0   dist_tumor < 0 
            figure(400)
            subplot(4,2,1+(i-1)*2)
            scatter(log2(double(cum_solid(ind_sub,2))),cum_solid(ind_sub,3),10,[0.8 0.8 0.8],'filled'); 
            hold on
            if sum(ind_sub) > 15 
                [f,x] = ksdensity([log2(double(cum_solid(ind_sub,2))) cum_solid(ind_sub,3)]);%,'PlotFcn','contour');
                [M,c] = contour(reshape(x(:,1),[30 30]),reshape(x(:,2),[30 30]),reshape(f,[30 30]),20);
                colormap('jet')
                c.LineWidth = 2;
            end
            ylim([0.2 1])
            xlim([2 8])
            title([title_for_plot ' Inside Tumor'])
            xlabel('log2 cell num')
% 
%             ind_sub = ind  & dist_tumor >= 0 & cum_solid(:,2) >= jumps(j) & cum_solid(:,2) < jumps(j+1) ; %  cum_ki67 cum_pd1(:,2) >= 0 dist_tumor >= 0
%             figure(400)
%             subplot(4,2,2+(i-1)*2)
%             scatter(log2(double(cum_solid(ind_sub,2))),cum_solid(ind_sub,3),10,[0.8 0.8 0.8],'filled'); 
%             hold on
%             if sum(ind_sub) > 15 
%                 [f,x] = ksdensity([log2(double(cum_solid(ind_sub,2))) cum_solid(ind_sub,3)]);%,'PlotFcn','contour');
%                 [M,c] = contour(reshape(x(:,1),[30 30]),reshape(x(:,2),[30 30]),reshape(f,[30 30]),20);
%                 colormap('jet')
%                 c.LineWidth = 2;
%             end
%             ylim([0.2 1])
%             xlim([2 8])
%             title([title_for_plot ' Outside Tumor'])
%             xlabel('log2 cell num')

    end
    
    
end

