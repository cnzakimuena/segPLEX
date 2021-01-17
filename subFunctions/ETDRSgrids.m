
function ETDRSgrids(patientsFolder)

%ETDRSGRIDS generates 2D and 3D ETDRS grids

% addpath(folderName1,...,folderNameN) adds the specified folders to the 
% top of the search path for the current MATLAB session
% genpath(folderName) returns a character vector containing a path name 
% that includes folderName and multiple levels of subfolders below 
% folderName (relative directory indicated by the dot symbol is collapsed 
% with the the specified folder, i.e., the subfolder 'subfunction')

%% list names of folders inside the patients folder

myDir = dir(patientsFolder);
dirFlags = [myDir.isdir] & ~strcmp({myDir.name},'.') & ~strcmp({myDir.name},'..');
nameFolds = myDir(dirFlags);

%% for each 3x3 subfolder, turn segmented data into network graph

for i = 1:numel(nameFolds)
    
    % assemble patient folder string
    folder = fullfile(patientsFolder, nameFolds(i).name);
    
    % add line to LOG
    disp(logit(folder, ['Initiating netGraph; ' nameFolds(i).name ' folder']))
    
    patientDir = dir(fullfile(folder, 'Results'));
    dirFlags = [patientDir.isdir] & ~strcmp({patientDir.name},'.') & ~strcmp({patientDir.name},'..');
    subFolders = patientDir(dirFlags);
    
    for k = 1:numel(subFolders)
        
        nameFold = subFolders(k).name;
        scanType = nameFold(1:2);
        eyeSide = nameFold(end-1:end);
        if strcmp(scanType, '3m')
            
            load(fullfile(folder,'Results', nameFold, 'segmentation.mat'))
            
            % visualise segmentation results
            %1- Without the segmentation
            %   figure;imshow3D(volumeStruc,[])
            %    %OR, for flow images
            %   figure;imshow3D(volumeFlow,[])
            %2- With the segmentation
            %   figure;imshow3D(volumeStruc,[],'plot',cat(3,RPEt,RPEb, RVIf, lBM),'LineWidth',2)
            %   %OR, for flow images
            %   figure;imshow3D(volumeFlow,[],'plot',cat(3,RPEt,RPEb, RVIf, lBM),'LineWidth',2)
            %3- If you want to see the volume in the nasal temporal direction (90degre
            %   rotation in the anteriot posterior direction)
            %   figure;imshow3D(permute(volumeStruc,[1,3,2]),[],'plot',permute(cat(3,RPEt,RPEb, RVIf, lBM),[2,1,3]),'LineWidth',2)
            
            % volume resize
            % volume size reduction factor for RAM limitation
            sizeRed = 600/1536;
            
            scanTag{1} = eyeSide;
            scanTag{2} = sizeRed;
            
            % Create 2D and 3D ETDRS grid directory
            if ~exist(fullfile(folder,'Results', nameFold, 'ETDRS_grid'), 'dir')
                mkdir(fullfile(folder,'Results', nameFold, 'ETDRS_grid'));
            end
            gridFolder = fullfile(folder,'Results', nameFold, 'ETDRS_grid');
            
            disp('begin dimAdjustAll')
            [vol_flow,vol_struc,BM,RVI] = dimAdjustAll(volumeFlow,volumeStruc,lBM,RVIf,sizeRed);
            disp('end dimAdjustAll')
            
            % volumeFlow orientation change to en-face direction
            enFace_Flow = zeros(size(vol_flow, 2), size(vol_flow, 3), size(vol_flow, 1));
            for ff = 1:size(vol_flow,1)
                disp(num2str(ff))
                enFace_im2 = mat2gray(reshape(vol_flow(ff,:,:), [size(vol_flow, 2), size(vol_flow, 3)]));
                enFace_Flow(:,:,ff) = enFace_im2;
            end
            % clearvars vol_flow
            
            % check data orientation integrity
            maxFlow = zeros(size(enFace_Flow, 1),size(enFace_Flow, 2));
            for h = 1:size(enFace_Flow, 1)
                disp(num2str(h))
                for hh = 1:size(enFace_Flow, 2)
                    maxFlow(h,hh) = max(enFace_Flow(h,hh,:));
                end
            end
            % clearvars enFace_Flow
            % rotation for consistency with fundus image orientation
            maxFlow = imrotate(maxFlow,-90);
            %                 currDimB = size(vol_struc, 2);
            %                 newDimV = 1536; % from 300
            %                 maxFlow =  imresize(maxFlow, [newDimV currDimB]);
            % figure;imshow(maxFlow,[])
            
            retinaMap = BM-RVI;
            % horizontal flip for consistency with fundus image orientation
            retinaMap = flip(retinaMap, 2);
            % figure;imtool(retinaMap,[])
            % figure;imagesc(retinaMap)
            
            % create the fundus size image
            [imageSizeY, imageSizeX] = size(retinaMap);
            [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
            % create the circle in the image
            centerIm = round([imageSizeX imageSizeY]/2); %image center
            centerX = centerIm(1);
            centerY = centerIm(2);
            % figure;imshow(retinaMap,[])
            % hold on
            % plot(centerX,centerY,'*r')
            % *foveal area mask (1.5x1.5mm)*
            diaDim = 1500; % real desired radius dimension, um
            radiusFac = 1.95/sizeRed; % conversion factor, um/px
            radius = round((diaDim/radiusFac)/2);
            circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
            retinaMapInner = retinaMap.*circlePixels;
            %figure; imshow(circlePixels,[])
            %figure; imshow(retinaMapInner,[])
            %figure; imshow([retinaMap retinaMapInner],[])
            
            [rowsOfMask, colsOfMask] = find(circlePixels); % foveal coordinates
            indOfMask = find(circlePixels); % foveal indexes
            minValue = min(retinaMap(circlePixels));
            indOfMin =  find(retinaMap(circlePixels) == minValue);
            rowsOfMin = rowsOfMask(indOfMin);
            colsOfMin = colsOfMask(indOfMin);
            
            % foveal center point
            fovCenterX = round(median(colsOfMin));
            fovCenterY = round(median(rowsOfMin));
            % figure;imshow(retinaMap,[])
            % figure;imshow(maxFlow,[])
            % hold on
            % plot(colsOfMin,rowsOfMin,'oc')
            % plot(fovCenterX,fovCenterY,'.r')
            
            % 'fundGrid' uses 'maxStruc' dimensions and foveal center coordinates to
            % generate ETDRS regions masks; 'regionsETDRS' 3D array has a volume
            % direction structure corresponding to:
            % [regionOD1, regionOD2, regionOD3, regionOD4, regionOD5, regionOD6, ...
            % regionOD7, regionOD8, regionOD9, grid]
            disp('begin fundRegions')
            regionsETDRS = fundRegions(maxFlow, fovCenterX, fovCenterY, sizeRed);
            disp('end fundRegions')
            %figure;imshow3D(regionsETDRS,[])
            % maxFlow(logical(regionsETDRS(:,:,1))) = 1;
            % figure;imshow(maxFlow,[])
            % hold on
            % plot(colsOfMin,rowsOfMin,'ob')
            % plot(fovCenterX,fovCenterY,'.r')
            
            %                 disp('begin fundGrid')
            %                 [gridETDRS, gridIm] = fundGrid(fundIm2D);
            %                 disp('end fundGrid')
            %                 %gridIm = fundIm2D(gridETDRS) = 255;
            %                 %figure; imshow(gridETDRS,[])
            %                 %figure; imshow(gridIm,[])
            
            maskZeros = zeros(size(vol_flow));
            for  v = 1:size(regionsETDRS,3)
                currentRegion = regionsETDRS(:,:,v);
                currVolume = maskZeros;
                for vv = 1:size(vol_flow,3)
                    currVolume(:,:,vv) = currentRegion;
                end
                structETDRS.regionsETDRS_3D{:,v} = logical(currVolume);
            end
            
            save(fullfile(folder,'Results', nameFold,'scanInfo.mat'),'scanTag');
            save(fullfile(folder,'Results', nameFold,'ETDRS_grid','2DregionsETDRS.mat'),'regionsETDRS','fovCenterX','fovCenterY');
            save(fullfile(folder,'Results', nameFold,'ETDRS_grid','3DregionsETDRS.mat'),'structETDRS');
            
        end
    end
end

                
                
                
                
                
                
                