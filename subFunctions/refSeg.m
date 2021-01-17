
function refSeg(patientsFolder)

%REFSEG performs refined segmentation prior to network graph implementation
%folderPath : path towards patient folders.
%All the data is saved in each patient folder in the subfolder 'Results'

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

%% for each 3x3 subfolder, segment converted PLEX data

for i = 1:numel(nameFolds) 
    
    % assemble patient folder string
    folder = fullfile(patientsFolder, nameFolds(i).name);
    
    try
        
        % add line to LOG
        disp(logit(folder, ['Initiating refSeg; ' nameFolds(i).name ' folder']))
        
        patientDir = dir(fullfile(folder, 'DataFiles'));
        dirFlags = [patientDir.isdir] & ~strcmp({patientDir.name},'.') & ~strcmp({patientDir.name},'..');
        subFolders = patientDir(dirFlags);
        
        for k = 1:numel(subFolders)
            
            nameFold = subFolders(k).name;
            scanType = nameFold(1:2);
            if strcmp(scanType, '3m')
                
                load(fullfile(folder,'DataFiles', nameFold, ['volume_aligned_' nameFold '.mat']))
                load(fullfile(folder,'DataFiles', nameFold, ['avgScans_' nameFold '.mat']))
                load(fullfile(folder,'DataFiles', nameFold, ['RegisteredImagesFlow_' nameFold '.mat']));
                load(fullfile(folder,'DataFiles', nameFold, ['RegisteredImages_' nameFold '.mat']));
                load(fullfile(folder,'DataFiles', nameFold, ['ImageList_' nameFold '.mat']), 'ImBrevis');
                
                % refine RPE segmentation (3D graph)
                disp('begin segmRPE ')
                [RPE1,RPE2] = segmRPE(avgScans, DeltaY);
                %RPE1 and RPE2 are respectively the coordinate of anterior and posterior RPE.
                %to see the results execute the following commented line:
                %figure;imshow3D(avgScans,[],'plot',cat(3,RPE1,RPE2),'LineWidth',2)
                disp('end segmRPE ')
                
                 % refine RVI segmentation (3D graph)
                disp('begin segmRVI ')
                [RVI] = segmRVI(volume_aligned,lRVIf, DeltaY);
                %RVI is a 2D map containing the coordinates of the RVI.
                %to see the results execute the following commented line:
                %figure;imshow3D(avgScans,[],'plot',cat(3,RVI),'LineWidth',2)
                disp('end segmRVI ')
                clear volume_aligned;
                
                %%% segment CSI (2* 2D graph) *NOT OPTIMIZED FOR PLEX Elite 9000*
                %disp('begin segmCSIv2 ')
                %[CSI]=segmentCSIv2({pathList},modelRf,'avgScans',avgScans,'ImageList',ImageList);
                %disp('end segmCSIv2 ')
                %%figure;imshow3D(avgScans,[],'plot',cat(3,CSI),'LineWidth',2)
                
                % unflattening for original volume and angiography segmentation
                
                disp('begin unflattening')
                [volumeStruc, RPEt,RPEb,RVIf] = unflattenVol(avgScans, RPE1, RPE2, RVI, err, trans);
                disp('end unflattening')
                volumeFlow = bscanstore2vol(bscanstoreFlow);
                %figure;imshow3D(volumeStruc,[],'plot',cat(3,RPEt,RPEb, RVIf, lBM),'LineWidth',2)
                %figure;imshow3D(volumeFlow,[],'plot',cat(3,RPEt,RPEb, RVIf, lBM),'LineWidth',2)
                
                if ~exist(fullfile(folder,'Results'), 'dir')
                    mkdir(fullfile(folder,'Results'));
                end
                if ~exist(fullfile(folder,'Results', nameFold), 'dir')
                    mkdir(fullfile(folder,'Results', nameFold));
                end
                
                save(fullfile(folder,'Results', nameFold, 'segmentationAlligned.mat'),'avgScans','RPE1','RPE2','RVI');
                save(fullfile(folder,'Results', nameFold, 'segmentation.mat'),'volumeStruc','volumeFlow','RPEt','RPEb','RVIf','lBM','-v7.3');
                
            end
        end
        
    catch exception
        errorString = ['Error in refSeg. Message:' exception.message buildCallStack(exception)];
        if ~exist(fullfile(pwd,'error'), 'dir')
            mkdir(fullfile(pwd,'error'));
        end
        disp(logit(fullfile(pwd, 'error'),errorString));
        continue
    end
    
end

disp(logit(folder,'Done refSeg'))        
        
        
        
        
        
        
        
        
        
        
        
        
        
        