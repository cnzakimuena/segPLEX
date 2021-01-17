
function preSeg(patientsFolder)

%PRESEG performs initial segmentation prior to refined segmentation
%folderPath : path towards patient folders.
%All the data is saved in each patient folder in the subfolder 'DataFiles'

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

%% for each 3x3 subfolder, pre-segment converted PLEX data

for i = 1:numel(nameFolds) 
    
    % assemble patient folder string
    folder = fullfile(patientsFolder, nameFolds(i).name);
    
    try

        % add line to LOG
        disp(logit(folder, ['Initiating preSeg; ' nameFolds(i).name ' folder']))
        
        patientDir = dir(fullfile(folder, 'ProcessedImages'));
        dirFlags = [patientDir.isdir] & ~strcmp({patientDir.name},'.') & ~strcmp({patientDir.name},'..');
        subFolders = patientDir(dirFlags);
        
        for k = 1:numel(subFolders)
            
            nameFold = subFolders(k).name;
            scanType = nameFold(1:2);
            if strcmp(scanType, '3m')

                [bscanstore, bscanstoreFlow] = bScanCrop(folder, nameFold); %crop the patient bscan (crop 5 colums on the left and on the right)
                load(fullfile(folder,'DataFiles', nameFold, ['RegisteredImagesFlow_' nameFold '.mat']));
                load(fullfile(folder,'DataFiles', nameFold, ['RegisteredImages_' nameFold '.mat']));
                load(fullfile(folder,'DataFiles', nameFold, ['ImageList_' nameFold '.mat']), 'ImBrevis')
                
                % create the volume : align all the bscan; flattening of the Bruch Membrane
                disp('begin bscanstore2volume3')

                [volume_aligned,lRVIf,lBM,lRPEt,lRPEb,err,trans] = bscanstore2volume3(bscanstore);
                disp('end bscanstore2volume3')
                clear bscanstore;
                
                % average bscans with adjacent bscans
                disp('begin lateralAverage')
                [avgScans,DeltaY]  = lateralAverage(volume_aligned, ImBrevis);
                disp('end lateralAverage')
                
                % conversion to uint8 to save memory space
                volume_aligned = uint8(mat2gray(volume_aligned)*255);
                avgScans = uint8(mat2gray(avgScans)*255);
                save(fullfile(folder,'DataFiles', nameFold, ['volume_aligned_' nameFold '.mat']),'volume_aligned','lRVIf','lBM','lRPEt','lRPEb','trans','err','-v7.3')
                save(fullfile(folder,'DataFiles', nameFold, ['avgScans_' nameFold '.mat']),'avgScans','DeltaY','-v7.3');
                % figure; imshow3D(volume_aligned,[])
                % figure; imshow3D(avgScans,[],'plot',cat(3,lRVIf,lRPEt,lRPEb))
                
            end
        end
        
    catch exception
        errorString = ['Error in preSeg. Message:' exception.message buildCallStack(exception)];
        if ~exist(fullfile(pwd,'error'), 'dir')
            mkdir(fullfile(pwd,'error'));
        end
        disp(logit(fullfile(pwd, 'error'),errorString));
        continue 
    end
    
end

disp(logit(folder,'Done preSeg'))
