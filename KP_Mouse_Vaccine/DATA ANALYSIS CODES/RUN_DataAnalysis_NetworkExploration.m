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
load([analfolder resufolder 'Results_CellType_' date '.mat'])
load([analfolder resufolder 'Results_Settings_' date '.mat'])
options.figOpt = 0;

filename.analfolder = analfolder;
filename.basefolder = basefolder;
filename.resufolder = resufolder;
options.date = date;

save([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'],'filename','options')
%%
LymphoNets.thresh.dist = 25; % 50
LymphoNets.thresh.netsize = 5; % 5;
cell_types = [11]; %  102 run this for both epithelial cells, lymphocytes and myeloid compartment

LymphoNets.Indexes   = uint16(zeros(size(MorpResults.Indexes)));
LymphoNets.Size      = uint16(zeros(size(MorpResults.Indexes)));
LymphoNets.DegreCent = uint16(zeros(size(MorpResults.Indexes)));
LymphoNets.ClnesCent = double(zeros(size(MorpResults.Indexes)));
LymphoNets.NetworkID = uint16(zeros(size(MorpResults.Indexes)));
LymphoNets.AllEdges  = uint16(zeros(size(MorpResults.Indexes)));
LymphoNets.LymphoEdges = uint16(zeros(size(MorpResults.Indexes)));

% loop around each tissue
for j = 1:length(options.MouseNum)    
    % find mouse number for plots only
    mousenum = options.MouseNum(j)
    
    % initialize variable for summary stats
    LymphoNets.Summary{j} = []; 
    
    % find cells from single tissue and their x-y coordinates
    index = MorpResults.Indexes == j & ROIResults.SpleenIndex == 0;
    coords = double([MorpResults.X(index) MorpResults.Y(index)]);
    
    
    % perform delaunay graph for all the cells
    S = delaunaygraph(coords,LymphoNets.thresh.dist,0);
    
    % calculate all the edges for each node for the full graph
    LymphoNets.AllEdges(index) = [uint16(full(sum(S))'); zeros(sum(index)-length(full(sum(S))),1)];
   
    
    % get the cell type desired ie lymphocytes & filter the nodes/edges 
    celltypefilt = CellType.Matrix(index,:);
    layer = ceil(cell_types(1)/10);
    index_t = celltypefilt(:,layer) == cell_types(1);
   
    S(~index_t,:) = 0;
    S(:,~index_t) = 0;
    T = graph(S);
    [bins,binsizes]=conncomp(T);
    nets = find(binsizes>LymphoNets.thresh.netsize);
    
    % calculate all the edges for each node for the full graph
    LymphoNets.LymphoEdges(index) = [uint16(full(sum(S))') ; zeros(sum(index)-length(full(sum(S))),1)];
    
    % initialize variables
    degree_cent = zeros(size(bins));
    closen_cent = zeros(size(bins));
    netsize =  zeros(size(bins));
    netId = zeros(size(bins));

    
    for n = 1:length(nets)
        % find nodes of subgraph
        S_sub = subgraph(T,bins==nets(n));
        
        
        
        % calculate centra
        degree_cent(bins==nets(n)) = centrality(S_sub,'degree');
        closen_cent(bins==nets(n)) = centrality(S_sub,'closeness');
        LymphoNets.Summary{j}.DegreCent(n) = median(centrality(S_sub,'degree'));
        LymphoNets.Summary{j}.ClnesCent(n) = median(centrality(S_sub,'closeness'));
        
        % calculate the size and composition
        netsize(bins==nets(n)) = binsizes(nets(n));
        netId(bins==nets(n)) = n;
        x_temp = coords(bins==nets(n),1);
        y_temp = coords(bins==nets(n),2);
        
        LymphoNets.Summary{j}.NumCells(n) = binsizes(nets(n));
        LymphoNets.Summary{j}.X_mean(n) = mean(x_temp);
        LymphoNets.Summary{j}.Y_mean(n) = mean(y_temp);
        LymphoNets.Summary{j}.Tcells(n) = sum(celltypefilt(bins==nets(n),3) == 111);
        LymphoNets.Summary{j}.Bcells(n) = sum(celltypefilt(bins==nets(n),3) == 112);
        LymphoNets.Summary{j}.Tregs(n) = sum(celltypefilt(bins==nets(n),4) == 1111);
        LymphoNets.Summary{j}.Thelps(n) = sum(celltypefilt(bins==nets(n),4) == 1112);
        LymphoNets.Summary{j}.Tcytox(n) = sum(celltypefilt(bins==nets(n),4) == 1113);
        
    end
    
    LymphoNets.Size(index) = padarray(netsize,[0 sum(index)-length(netsize)],'post');
    LymphoNets.DegreCent(index) = padarray(degree_cent,[0 sum(index)-length(degree_cent)],'post');
    LymphoNets.ClnesCent(index) = padarray(closen_cent,[0 sum(index)-length(closen_cent)],'post');
    LymphoNets.NetworkID(index) = padarray(netId,[0 sum(index)-length(netId)],'post');
end

save([filename.analfolder filename.resufolder 'Results_Nets_' options.date '_dist' num2str(LymphoNets.thresh.dist) '.mat'],'LymphoNets')

% %%
% clear all
% codedir = 'Y:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\CycIF Codes\Utility Functions';
% addpath(codedir)
% filename.basefolder = 'Y:\sorger\data\IN_Cell_Analyzer_6000\Claire\2020_01_IFNgR_Mouse_Lung\';
% filename.analfolder = [filename.basefolder 'ANALYSIS\'];
% filename.resufolder = 'Analysis_Results_20200709\'; 
% options.date = '20201024';
% 
% load([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'])
% load([filename.analfolder filename.resufolder 'Results_ROI_' options.date '.mat'])
% % load([filename.analfolder filename.resufolder 'Results_Norm_' options.date '.mat']);
% % load([filename.analfolder filename.resufolder 'Results_Filt_' options.date '.mat']);
% load([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'])
% load([filename.analfolder filename.resufolder 'Results_CellType_' options.date '.mat'])
% % load([filename.analfolder filename.resufolder 'Results_Nets_' options.date '.mat'])
% load([filename.analfolder filename.resufolder 'Results_Nets_' options.date '_dist25.mat'])

% %% calculate summary metric for the nets
% 
% allnets.size  = [];
% allnets.bfrac = [];
% allnets.tfrac = [];
% allnets.tRfrac = [];
% allnets.tHfrac = [];
% allnets.tCfrac = [];
% allnets.mouse = [];
% allnets.mousegroup = [];
% for j = 1:length(options.MouseNum)
%     allnets.size = [allnets.size; LymphoNets.Summary{j}.NumCells'];
%     allnets.bfrac = [allnets.bfrac; LymphoNets.Summary{j}.Bcells'./LymphoNets.Summary{j}.NumCells'];
%     allnets.tfrac = [allnets.tfrac; LymphoNets.Summary{j}.Tcells'./LymphoNets.Summary{j}.NumCells'];
%     allnets.tRfrac = [allnets.tRfrac; LymphoNets.Summary{j}.Tregs'./LymphoNets.Summary{j}.NumCells'];
%     allnets.tHfrac = [allnets.tHfrac; LymphoNets.Summary{j}.Thelps'./LymphoNets.Summary{j}.NumCells'];
%     allnets.tCfrac = [allnets.tCfrac; LymphoNets.Summary{j}.Tcytox'./LymphoNets.Summary{j}.NumCells'];
%     allnets.mouse = [allnets.mouse; zeros(length(LymphoNets.Summary{j}.NumCells),1)+options.MouseNum(j)];
%     allnets.mousegroup = [allnets.mousegroup; zeros(length(LymphoNets.Summary{j}.NumCells),1)+options.MouseGroup(j)];
%     
% end
% 
% %%
% colorbygroup = {'b','c','r','m','g'};
% cum_nets = cell(1,4);
% for j = 1:28
%     index = MorpResults.Indexes == j;
%     if options.MouseGroup(j) ==  3
%         
%         figure(5000+j)
%         scatter(MorpResults.X(index & ROIResults.TumorIndex > 0),MorpResults.Y(index  & ROIResults.TumorIndex > 0 ),10,[255 153 51]/255,'filled')
%         set(gca,'Ydir','reverse')
%         hold on
%         scatter(MorpResults.X(index & ROIResults.TumorIndex == 0),MorpResults.Y(index  & ROIResults.TumorIndex == 0 ),10,[51 255 255]/255,'filled','MarkerFaceAlpha',1)%0.05)
%         scatter(MorpResults.X(index & LymphoNets.NetworkID > 0),MorpResults.Y(index  & LymphoNets.NetworkID > 0),10,[255 51 255]/255,'filled')
%         
%         viscircles([LymphoNets.Summary{j}.X_mean' LymphoNets.Summary{j}.Y_mean'],LymphoNets.Summary{j}.NumCells'+0.1,'Color','k','LineStyle','--');
%         viscircles([LymphoNets.Summary{j}.X_mean' LymphoNets.Summary{j}.Y_mean'],LymphoNets.Summary{j}.Tcells'+0.1,'Color','b');
%         viscircles([LymphoNets.Summary{j}.X_mean' LymphoNets.Summary{j}.Y_mean'],LymphoNets.Summary{j}.Bcells'+0.1,'Color','r');
%        
%     end
%     cum_nets{options.MouseGroup(j)} = [cum_nets{options.MouseGroup(j)} LymphoNets.Summary{j}.NumCells];
% end
% %%
% th_cell = 15;
% 
% for i = 1:4
%     mean(cum_nets{i}(cum_nets{i}>th_cell))
%     sum(cum_nets{i}>th_cell)
%     figure(11)
%     subplot(1,4,i)
%     [n,h]=hist(cum_nets{i}(cum_nets{i}>th_cell),linspace(th_cell,th_cell*10,10));
%     bar(h,n)%/sum(n))
%     hold on
%     ylim([0 max([80 n])*1.02])
% end
% suptitle('Distribution of nets by size')
% 
% %%% COMPOSITION
% mean_Bfrac = [];
% mean_Tfrac = [];
% mean_tCfrac = [];
% mean_tHfrac = [];
% mean_tRfrac = [];
% tot= [];
% x = [];
% ylim_cum = 0.7;
% 
% bins = 40;
% jump = 10;
% max_net = bins*jump;
% % find
% for j = 5
%     for i = 1:bins
%         % find nets of certain size
%         ind_size = allnets.size > jump*(i-1) & allnets.size <= jump*i & allnets.mousegroup < j;
%         mean_Bfrac(i,:) = [nanmean(allnets.bfrac(ind_size)) prctile(allnets.bfrac(ind_size),25) prctile(allnets.bfrac(ind_size),75) ]; 
%         mean_Tfrac(i,:) = [nanmean(allnets.tfrac(ind_size)) prctile(allnets.tfrac(ind_size),25) prctile(allnets.tfrac(ind_size),75) ]; 
%         mean_tCfrac(i,:) = [nanmean(allnets.tCfrac(ind_size)./allnets.tfrac(ind_size)) prctile(allnets.tCfrac(ind_size)./allnets.tfrac(ind_size),25) prctile(allnets.tCfrac(ind_size)./allnets.tfrac(ind_size),75) ]; 
%         mean_tHfrac(i,:) = [nanmean(allnets.tHfrac(ind_size)./allnets.tfrac(ind_size)) prctile(allnets.tHfrac(ind_size)./allnets.tfrac(ind_size),25) prctile(allnets.tHfrac(ind_size)./allnets.tfrac(ind_size),75) ]; 
%         mean_tRfrac(i,:) = [nanmean(allnets.tRfrac(ind_size)./allnets.tfrac(ind_size)) prctile(allnets.tRfrac(ind_size)./allnets.tfrac(ind_size),25) prctile(allnets.tRfrac(ind_size)./allnets.tfrac(ind_size),75) ];  
%         x(i,1) = i*jump - jump/2;
%         tot(i) = nansum(ind_size);
%         if tot(i) < 10
%             mean_Bfrac(i,:) = [ NaN NaN NaN];
%             mean_Tfrac(i,:) = [ NaN NaN NaN];
%             mean_tCfrac(i,:) = [ NaN NaN NaN];
%             mean_tHfrac(i,:) = [ NaN NaN NaN];
%             mean_tRfrac(i,:) = [ NaN NaN NaN];
%         end
%     end
%     figure(2021)
%     subplot(1,2,1)
%     scatter(x-0.5,mean_Bfrac(:,1),log2(tot+1.1)*3,'r','filled')
%     hold on
%     scatter(x+0.5,mean_Tfrac(:,1),log2(tot+1.1)*3,'b','filled')
%     plot([x-.5 x-.5]',mean_Bfrac(:,2:3)','-r')
%     
%     plot([x+.5 x+.5]',mean_Tfrac(:,2:3)','-b')
%     legend('B cell frac','T cell frac')
%     ylim_cum = max([ylim_cum max([mean_Bfrac,mean_Tfrac])]);
%     ylim([0 ylim_cum*1.05])
%     
%     subplot(1,2,2)
%     scatter(x,mean_tCfrac(:,1),5*3,'filled')
%     hold on
%     scatter(x,mean_tHfrac(:,1),5*3,'filled')
%     scatter(x,mean_tRfrac(:,1),5*3,'filled')
%     legend('Tc cell frac','Th cell frac','Tr cell frac')
%     ylim([0 ylim_cum*1.05])
% end
% 
% %%
% forplot.grouplist = {'Cre+dgTom','Cre+dgCxcl10','LucOS+dgTom','LucOS+dgCxcl10'};
% forplot.group_mat = [266 267 268 269 270; 
%              271 272 273 274 275; 
%              276 277 278 279 280; 
%              281 282 283 284 285];
%          
% bin_size  = 1000;
% min_cell = 50;
% 
% % loop around each mouse
% for g = 1:3%size(forplot.group_mat,1)
%     B_N = [];
%     % loop around the mouse within each condition
%     for i = 1:size(forplot.group_mat,2)
%         % isolate cells of specific mouse
%         rois_ind = find(options.MouseNum==forplot.group_mat(g,i));
%         
%         for r = 1:length(rois_ind)
%             index = MorpResults.Indexes == rois_ind(r) & ROIResults.SpleenIndex == 0;
%             cent_x = double(ceil(MorpResults.X(index)/bin_size));
%             cent_y = double(ceil(MorpResults.Y(index)/bin_size));
%             
%             edges{1} = linspace(0,max(cent_x),max(cent_x)+1)+0.5;
%             edges{2} = linspace(0,max(cent_y),max(cent_y)+1)+0.5;
%             [N,c] = hist3([cent_x cent_y],'Edges',edges);
%             
%             index = MorpResults.Indexes == rois_ind(r) & ROIResults.SpleenIndex == 0 & CellType.Matrix(:,3) == 111 ;% & LymphoNets.Size <  5;
%             cent_x = double(ceil(MorpResults.X(index)/bin_size));
%             cent_y = double(ceil(MorpResults.Y(index)/bin_size));
%             [B,c] = hist3([cent_x cent_y],'Edges',edges);
%             
%             B_N_temp = B(N(:)>min_cell)./N(N(:)>min_cell);
%             B_N = [B_N; B_N_temp];
%           
%         end
%         
% %         
% %         index = ismember(MorpResults.Indexes,rois_ind) & ROIResults.SpleenIndex == 0;
% %         
% %         bcells_tot(g,i) = sum(CellType.Matrix(index,3) == 112) / sum(CellType.Matrix(index,1) == 1); 
% %         
% %         bcells_free(g,i) = sum(CellType.Matrix(index,3) == 112  & LymphoNets.Size(index) == 0)  / sum(CellType.Matrix(index,1) == 1); 
% %         
% %         index = ismember(MorpResults.Indexes,rois_ind) & ROIResults.SpleenIndex == 0 & ROIResults.TumorIndex == 0;
% %       
% %         bcells_free(g,i) = sum(CellType.Matrix(index,3) == 112 & LymphoNets.Size(index) < 15)  ;%/ sum(CellType.Matrix(index,2) == 11);
%         
%     end
%     figure(2)
%     [n,h] = hist(B_N,linspace(0,0.3,61));
%     plot(h,n/sum(n))
%     xlim([0 .3])
%     hold on
% end
% %%
% plotvect = [ones(5,1)-(6-(1:5)'-3)*0.1; 2*ones(5,1)-(6-(1:5)'-3)*0.1; 3*ones(5,1)-(6-(1:5)'-3)*0.1; 4*ones(5,1)-(6-(1:5)'-3)*0.1];
% 
% temp_mat = bcells_tot';
% median_forbar = median(temp_mat);
% figure(1)
% subplot(1,2,1)
% bar([1:4],median_forbar,0.2,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
% hold on
% scatter(plotvect,temp_mat(:),40,round(plotvect),'filled')
% caxis([0.5 4])
% xlim([0.5 4.5])
% ylim([0 max(temp_mat(:))*1.1])
% colormap('jet')
% set(gca,'XTick',1:4,'XTickLabel',forplot.grouplist,'XTickLabelRotation',45)
% title('B cells per mouse')
% ylabel('% of all cells')
% 
% temp_mat = bcells_free';
% median_forbar = median(temp_mat);
% figure(1)
% subplot(1,2,2)
% bar([1:4],median_forbar,0.2,'EdgeColor','k','FaceColor',[0.9 0.9 0.9])
% hold on
% scatter(plotvect,temp_mat(:),40,round(plotvect),'filled')
% caxis([0.5 4])
% xlim([0.5 4.5])
% ylim([0 max(temp_mat(:))*1.1])
% colormap('jet')
% set(gca,'XTick',1:4,'XTickLabel',forplot.grouplist,'XTickLabelRotation',45)
% title('free B cells per mouse')
% ylabel('% of all cells')
% 
% %%
% 
% tsum = allnets.tRfrac + allnets.tHfrac + allnets.tCfrac;
% 
% %
% th = 25;
% figure;
% for j = [1 3]
%     subplot(1,5,1)
%     [n,h]=hist(allnets.bfrac(allnets.mousegroup == j & allnets.size > th),linspace(0,1,11));
%     plot(h,n/sum(n))
%     hold on
%     subplot(1,5,2)
%     [n,h]=hist(allnets.tfrac(allnets.mousegroup == j & allnets.size > th),linspace(0,1,11));
%     plot(h,n/sum(n))
%     hold on
%     subplot(1,5,3)
%     [n,h]=hist(allnets.tRfrac(allnets.mousegroup == j & allnets.size > th),linspace(0,1,11));
%     plot(h,n/sum(n))
%     hold on
%     subplot(1,5,4)
%     [n,h]=hist(allnets.tHfrac(allnets.mousegroup == j & allnets.size > th),linspace(0,1,11));
%     plot(h,n/sum(n))
%     hold on
%     subplot(1,5,5)
%     [n,h]=hist(allnets.tCfrac(allnets.mousegroup == j & allnets.size > th),linspace(0,1,11));
%     plot(h,n/sum(n))
%     hold on
% end
% 
% 
% %%
% mean_Bfrac = [];
% tot= [];
% x = [];
% 
% for i = 1:max(allnets.size)
%     if sum(allnets.size==i) > 2
%     mean_Bfrac(i) = mean(allnets.bfrac(allnets.size==i));
%     x(i) = i;
%     tot(i) = sum(allnets.size==i);
%     end
% end
% figure(2020)
% scatter(x,mean_Bfrac,log2(tot+1.1)*2,'filled')
% % hold on
% % bounds = [5 60; 60 120];
% % for i = 1%:size(bounds,1)
% %     [p,S]=polyfit(allnets.size(allnets.size>=bounds(i,1) & allnets.size<bounds(i,2)),allnets.bfrac(allnets.size>=bounds(i,1) & allnets.size<bounds(i,2)),1);
% %     x0 = bounds(i,1):bounds(i,2);
% %     plot(x0,p(2)+p(1)*x0)
% % end
% ylabel('B cell fraction')
% xlabel('Number of cells in network')
% 
% 
% 







