function PreProcess_Step4_ExtractROIs(filename,options)

% Retrieve ROI's from ImageJ scoring
% basefolder = filename.basefolder; % 'Claire\PD-1_Checkpoint_Mouse_Lung_03-2020\';
% analfolder = filename.analfolder; % [ImStor basefolder 'ANALYSIS\'];
% roi_folder = filename.roifolder; % 'ROIs\';
% mont_folder = filename.montfolder; % 'MontageforROI_Lv3\';
% mont_suffix = filename.montsuffix ; % '_montage.tif';

% folders = filename.folders ; %{ 'MLB_1059','MLB_1060','MLB_1061','MLB_1062','MLB_1063','MLB_1064', ...
%             'MLB_1065','MLB_1066','MLB_1067','MLB_1068','MLB_1069','MLB_1070', ...
%             'MLB_1079','MLB_1080','MLB_1081','MLB_1082','MLB_1083','MLB_1084', ...
%             'MLB_1085','MLB_1086','MLB_1087','MLB_1088','MLB_1089','MLB_1090', ...
%                     };  
                
% montagelevel = options.ROI.pyramidlevel; % 3;
% 
% ROI_names = options.PathLabels; %{'Tumors','Spleen'};%
% 
% tumor_count = 0;
% TumorTrack = [];

for fol = 1:length(filename.folders)
    
    montfile = [filename.analfolder filename.montfolder filename.folders{fol} filename.montsuffix];

    montinfo = imfinfo(montfile); 
    mcolsize = montinfo.Width;
    mrowsize = montinfo.Height; 

    % first do the spleen
    % create ROI image
    
    for p = 1:length(options.PathLabels)
        ROI_img = zeros(mrowsize,mcolsize);
        try
            ROIfile = [filename.analfolder filename.roifolder filename.folders{fol} '_' options.PathLabels{p} '.roi'];
            ROI_mont = ReadImageJROI(ROIfile);
        catch
            try
                ROIfile = [filename.analfolder filename.roifolder filename.folders{fol} '_' options.PathLabels{p} '.zip'];
                ROI_mont = ReadImageJROI(ROIfile);
            catch
                disp([filename.analfolder filename.roifolder filename.folders{fol} '_' options.PathLabels{p} ' not found'])
                continue
            end
        end 
        
        for i3 = 1:length(ROI_mont)
            if length(ROI_mont) ==1
                coords = ROI_mont.mnCoordinates;
            else
                coords = ROI_mont{1,i3}.mnCoordinates;
            end
            coords = [coords(end,:); coords];
            ROI_img(sub2ind([mrowsize mcolsize],coords(:,2),coords(:,1))) = 1;
            %Connecting two points in the matrix and making connection = 1
            for j1 = 1:length(coords)-1
                xline = round(linspace(coords(j1,2), coords(j1+1,2),100));
                yline = round(linspace(coords(j1,1), coords(j1+1,1),100));
                idx=sub2ind([mrowsize mcolsize],xline,yline);
                ROI_img(idx)=i3;
            end
        end
        ROI_img = imfill(ROI_img,'holes');
        
        saveastiff(uint8(ROI_img),[filename.analfolder filename.roifolder filename.folders{fol} '_' options.PathLabels{p} '.tif'])
    
        figure(fol)
        subplot(1,length(options.PathLabels),p)
        imshow(uint8(ROI_img),[])
    end
    
%     
%     ROI_img = zeros(mrowsize,mcolsize);
%     
%     Spleen_ROIfile = [filename.analfolder filename.roifolder filename.folders{fol} '_Spleen.roi'];
%     Spleen_ROI_mont = ReadImageJROI(Spleen_ROIfile);
%     coords = Spleen_ROI_mont.mnCoordinates;
%     % check that there is no out of bounds coordinates
%     coords(coords==0) = 1;
%     coords = [coords(end,:); coords];
%     ROI_img(sub2ind([mrowsize mcolsize],coords(:,2),coords(:,1))) = 1;
%     for j1 = 1:length(coords)-1
%         xline = round(linspace(coords(j1,2), coords(j1+1,2),100));
%         yline = round(linspace(coords(j1,1), coords(j1+1,1),100));
%         idx=sub2ind([mrowsize mcolsize],xline,yline);
%         ROI_img(idx)=1;
%     end
%     
%     Spleen_Img = imfill(ROI_img,'holes');
%     % now dilate the spleen image to be sure to get all of it
%     temp = Spleen_Img;
%     Spleen_Img = imdilate(Spleen_Img,strel('disk',50));
%     Spleen_Img = imerode(Spleen_Img,strel('disk',50));
%     
%     % then do tumors
%     ROI_img = zeros(mrowsize,mcolsize);
%     
%     Tumors_ROIfile = [filename.analfolder filename.roifolder filename.folders{fol} '_Tumors.zip'];
%     Tumors_ROI_mont = ReadImageJROI(Tumors_ROIfile);
%     
%     for i3 = 1:length(Tumors_ROI_mont)
%         if length(Tumors_ROI_mont) ==1
%             coords = ROI_mont.mnCoordinates;
%         else
%             coords = Tumors_ROI_mont{1,i3}.mnCoordinates;
%         end
%         coords = [coords(end,:); coords];
%         ROI_img(sub2ind([mrowsize mcolsize],coords(:,2),coords(:,1))) = 1;
%         %Connecting two points in the matrix and making connection = 1
%         for j1 = 1:length(coords)-1
%             xline = round(linspace(coords(j1,2), coords(j1+1,2),100));
%             yline = round(linspace(coords(j1,1), coords(j1+1,1),100));
%             idx=sub2ind([mrowsize mcolsize],xline,yline);
%             ROI_img(idx)=1;
%         end
%     end
%     Tumor_Img = imfill(ROI_img,'holes');
%     
%     saveastiff(uint8(Spleen_Img),[analfolder roi_folder folders{fol}(5:8) '_Spleen.tif'])
%     saveastiff(uint8(Tumor_Img),[analfolder roi_folder folders{fol}(5:8) '_Tumors.tif'])
    
%     figure
%     subplot(1,2,1)
%     imshow(Tumor_Img,[])
%     subplot(1,2,2)
%     imshow(Spleen_Img,[])
end