
function [bscanstore, bscanstoreFlow] = bScanCrop(folderPath, foldID)

%BSCANCROP crop the original bscans (crop 5 colums on the left and 5 colums
%on the right. 
%   VARARGIN is the patient folder.
%   FOLDID is the the identifier corresponding to the bscan set.
%   BSCANSTORE is a cell containing all bscan cropped. 
%   (warning : as bscanstore are generated after the IMAGELIST table (see
%   convertPLEX) the spatial data in the IMAGELIST such as the
%   position of the bscans on the fundus image doesn't consider the
%   cropping). 

%dirlist = adaptToHMRpath(varargin{1});

disp(logit(folderPath,'Starting mapPseudoRegistration'))
cropNum = 5;%cropping limits specification

try

    load(fullfile(folderPath, 'DataFiles', foldID, ['ImageList_' foldID '.mat']), 'ImBrevis');    
    if exist(fullfile(folderPath, 'DataFiles', foldID, ['TrimInfo_' foldID '.txt']),'file')
        cropLimits = dlmread(fullfile(folderPath, 'DataFiles', foldID, ['TrimInfo_' foldID '.txt']),',');
    else
        fname = fullfile(folderPath,'DataFiles', foldID, ['TrimInfo_' foldID '.txt']);
        col   = [cropNum;ImBrevis.width(1)-cropNum];%cropping limits implementation
        dlmwrite(fname,min(max(1,sort(col)),ImBrevis.width(1))','precision','%g','delimiter',',')
        cropLimits = dlmread(fullfile(folderPath, 'DataFiles', foldID, ['TrimInfo_' foldID '.txt']),',');
    end    
    
    pngList1 = dir(fullfile(folderPath,'ProcessedImages', foldID, [foldID '_cube_z'], '*.png'));
    numframes1  = numel(pngList1);
    bscanstore = cell(numframes1,1);
    for frame1 = 1:numframes1
        bscan = imread(fullfile(folderPath,'ProcessedImages', foldID, [foldID '_cube_z'], pngList1(frame1).name));
        if exist('cropLimits','var')
            bscanstore{frame1} = bscan(:,cropLimits(1):cropLimits(2));
        else
            bscanstore{frame1} = bscan(:,cropNum:size(bscan,2)-cropNum);%cropping limits implementation
        end
    end
    %bscanAlign =alignBscan(bscanstore);
    skippedind = [];
    start      = 1;
    sizeStore = size(bscanstore{1});
    save(fullfile(folderPath,'DataFiles', foldID, ['RegisteredImages_' foldID '.mat']),'bscanstore','skippedind','start','sizeStore');
    %createAllDatawp(folder)
    
    pngList2 = dir(fullfile(folderPath,'ProcessedImages', foldID, [foldID '_FlowCube_z'], '*.png'));
    numframes2  = numel(pngList2);    
    bscanstoreFlow = cell(numframes2,1);
    for frame2 = 1:numframes2
        bscanFlow = imread(fullfile(folderPath,'ProcessedImages', foldID, [foldID '_FlowCube_z'], pngList2(frame2).name));
        if exist('cropLimits','var')
            bscanstoreFlow{frame2} = bscanFlow(:,cropLimits(1):cropLimits(2));
        else
            bscanstoreFlow{frame2} = bscanFlow(:,cropNum:size(bscanFlow,2)-cropNum);
        end
    end
    skippedindFlow = [];
    startFlow      = 1;
    sizeFlowStore = size(bscanstoreFlow{1});
    save(fullfile(folderPath,'DataFiles', foldID, ['RegisteredImagesFlow_' foldID '.mat']),'bscanstoreFlow','skippedindFlow','startFlow','sizeFlowStore');
 
catch exception
    
    errorString = ['Error in mapPseudoRegistration. In' folderPath ' Message:' exception.message buildCallStack(exception)];
    disp(logit(folderPath,errorString));
    
end
                             
disp(logit(folderPath,'Done mapPseudoRegistration'))
    
end



