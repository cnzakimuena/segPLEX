
function [volumeFlow,volumeStruc,BMf,RVIf] = dimAdjustAll(vol_flow,vol_struc,BM,RVI,sizeReducFac)

% modify 3D volume aspect ratio to true proportions
% dimension decreased by reduction factor
decNewDimA = round(size(vol_flow, 1)*sizeReducFac); 
decNewDimB = round(size(vol_flow, 2)*sizeReducFac);
% assuming same real distance in a-scan and volume direction (3mmx3mm) 
decNewDimV = decNewDimA; 
% returns the volume B that has the number of rows, columns, and planes 
% specified by the three-element vector [numrows numcols numplanes].
volumeFlow = imresize3(vol_flow,[decNewDimA decNewDimB decNewDimV]);
volumeStruc = imresize3(vol_struc,[decNewDimA decNewDimB decNewDimV]);
% figure;imshow3D(avgScans,[])

% modify 2D segmentation surfaces aspect ratio to true proportions
BMf =  round(imresize(BM, [decNewDimV decNewDimB]).*sizeReducFac);
RVIf = round(imresize(RVI, [decNewDimV decNewDimB]).*sizeReducFac);
% trans = round(imresize(trans, [newDimV newDimB]));

end


