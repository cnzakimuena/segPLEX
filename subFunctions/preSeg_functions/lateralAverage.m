function [avgScans,DeltaY] = lateralAverage(volume_aligned, ImageList)
%%LATERALAVERAGE do a weighted average for each bscan with its adjacent bscans in order to reduce the speckle noise.
%   VOLUME_ALIGNED : OCT volume given by the function bscanstore2volume3
%   IMAGELIST      : data for each bscan (given by the function convertPLEX)
%   AVGSCANS       : average volume
%   DELTAY         : distance between 2 bscan.
%   DeltaY is an important data because it tunes the weights of the means.
%   If 2 adjacent bscans are far, the weights of the mean are weak and reciprocally

start = 1;
volume_aligned = single(mat2gray(volume_aligned));

strCell = ImageList{start, 'protocol'};
refStr = strCell{1};
imgType = refStr(8:9);                                    
if strcmp(imgType, '3m')
    [bscanDim, volDim] = deal(3);
elseif strcmp(imgType, '6m')
    [bscanDim, volDim] = deal(6);
elseif strcmp(imgType, '9m')
    [bscanDim, volDim] = deal(9);
elseif strcmp(imgType, '12')
    [bscanDim, volDim] = deal(12);
elseif strcmp(imgType, '15')
    bscanDim = 9;
    volDim = 15;    
end
[nl,nc,numframes] = size(volume_aligned);
% b-scan lateral resolution in mm/pixel 
DeltaX = bscanDim / ImageList{start, 'width'}; 
% volume direction resolution from b-scan spacing in mm/pixel
DeltaY = volDim / numframes; 

p=5;
SigmaFilterScans = max(1,ceil(p * DeltaX / DeltaY)); %nb of bscan in the average from either side
interScansFilter = exp(-(-SigmaFilterScans:SigmaFilterScans).^2 / 2 / SigmaFilterScans^2); %weight of the gaussian
interScansFilter = interScansFilter / sum(interScansFilter);%normalisation
avgScans = single(zeros(nl,nc,numframes));

for frame = 1:numframes%indToProcess
    try
        
        startFrame = max(1,frame-SigmaFilterScans);
        lastFrame  = min(numframes,frame+SigmaFilterScans);
        
        allAux = [];

        % Concatenate images to average
        for avgFrame = startFrame:lastFrame
            avgWeight = interScansFilter(avgFrame - frame + SigmaFilterScans + 1);
            allAux = cat(3,allAux,volume_aligned(:,:,avgFrame) * avgWeight);
        end
        
        % Compute weighted average
        avgScans(:,:,frame) = nansum(allAux,3) / sum(interScansFilter((startFrame:lastFrame) - frame + SigmaFilterScans + 1));
        disp(frame)
    catch
        disp(logit(savedir,['Error preProcessFrames at frame:' num2str(frame)]));
    end
end
%figure; imshow3D(avgScans,[])
%figure; imshow3D(volume_aligned,[])

end