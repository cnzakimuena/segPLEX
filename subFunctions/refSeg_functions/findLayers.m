function layers = findLayers(labels,s1,s2,s3,s4)
%FINDLAYERS transforms the labels into a map of coordinates
%   s1,s2,s3,s4 are respectively the number of row, colums (A-scan), Bscans
%   and surfaces.
%   LABELS : label of each pixels in the 4d volume (given by
%   graphCutPierre2)
%   LAYERs is a cell containing  the coordinate map for each layer
%   segmented
imageResults = zeros(s1,s2,s3,s4,'single');
imageResults(find(labels)) = 1; %transform the labels vector into a 4D volume to obtain the segmentation in 3+1D
imageResults(1,:,:,:) = 1;
%figure; imshow3D(imageResultsRPE,[])
layers = cell(1,s4);
for i =1:s4 %calculate the gradient to see the interface
    layers_i = nan(s3,s2);
    imageResults_i = imfilter(imageResults(:,:,:,i),[1;-1]);
    ind_i = find(imageResults_i==1);
    [l,c,z] = ind2sub(size(imageResults_i),ind_i);
    layers_i(sub2ind([s3,s2],z,c)) = l;
    layers{i} = layers_i;
end
%figure; subplot(2,1,1),imshow(imageResults_i(:,:,1),[]);subplot(2,1,2),imshow(imageResults(:,:,1,i),[])
    


end