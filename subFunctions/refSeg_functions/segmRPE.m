function [RPE1,RPE2] = segmRPE(volume_aligned,DeltaY)
%SEGMRPE segments the top and bottom of the RPE with a 3D graph
%   volume_aligned : 3D OCT volume
%   DeltaY  : distance between 2 bscans
%   RPE1 : anterior RPE. 2D matrices where each lines are the ordinate of
%   the anterior RPE for each bscan
%   RPE2 : posterior RPE 2D matrices where each lines are the ordinate of
%   the posterior RPE for each bscan
%   for more details of the 3D graph algorithm,see \readme\Algorithme de segmentation des couches rétiniennes et choroïdiennes.docs) 

[~,~,nbz] = size(volume_aligned); %nb of bscan
ref =510; %position of the BM
wt = 60; wb = 20; 
volRPE =single(mat2gray(volume_aligned(ref-wt:ref+wb,:,:))); %we only select a volume around the RPE
volRPEf = imfilter(volRPE,fspecial('gaussian',[3,3],1),'symmetric'); %gaussian filter
%figure; imshow3D(volRPEf );

volRPEf2=imclose(volRPEf,strel('disk',8,8)); %closure : dilatation and erosion. to homogenize the RPE 
%figure; imshow3D(volRPEf2,[]);
[inflex1,inflex2,~,~] = inflexImgAngles(volRPEf2,[0]); %find the position of inflexion points and its weights (see \readme\Algorithme de segmentation des couches rétiniennes et choroïdiennes.docs) 

neg = inflex1==0; inflex1 = inflex1-10*neg; %impose a negative weight if there is no inflexion point
neg = inflex2==0; inflex2 = inflex2-10*neg;

%figure; imshow3D(inflex1 );figure; imshow3D(inflex2 );figure; imshow3D(inflex3,[] );
% figure; imshow3D(grad );

deltaMin=[2];deltaMax=[50]; %deltaMin :minimum distance between anterior and posterior RPE
                            %deltaMax :maximum distance between anterior and posterior RPE 
jumpMaxCol = 2;% the maximum jump a layer can do between each Ascan is 'jumpMaxCol' pixels 
if DeltaY >0.07
    jumpMaxBscan =5;%the maximum jump a layer can do between each Bscan line is 'jumpMaxBscan' pixels
                    %this number is dependant of the distance between each
                    %bscan (DeltaY). 
else
    jumpMaxBscan = 3;
end
%figure; imshow3D(inflex1,[])
%figure; imshow3D(inflex2,[])

%graph Cut
% if B-scans are larger than design (512x496), graph will be excuted in two
% parts, see '(mod)' below
if size(volume_aligned,2) > 492
    nbBscan = floor(nbz/10); % change dividing number according value of size(volume_aligned,2)
else
    nbBscan = nbz;
end

RPE1 = [];
RPE2 = [];

for i =1:floor(nbz/nbBscan)% this for loop is useless because there is only one loop.
    % (mod) You can modify the value of nbBscan if you want to divide the volume in slices if your RAM is not big enough
    
    disp(['begin slice ' num2str(i)])
    
    if length((i-1)*nbBscan+1:min((i+1)*nbBscan,nbz))<nbBscan*2
        inter = (i-1)*nbBscan+1:min((i+1)*nbBscan,nbz);
    else
        inter = (i-1)*nbBscan+1:min((i)*nbBscan,nbz);
    end
%the following variable voxel is a 4D matrice. In its function, it the
%concatenation of the weights for the anterior RPE (inflex1) and the
%weights for the posterior RPE (inflex2)
voxels  = cat(4,inflex1(:,:,inter)*255,inflex2(:,:,inter)*255);
voxels (isnan(voxels)) = 0;

%segmentation (3D graph): labels 
[~, labels] = graphCutPierre2(voxels,jumpMaxCol,jumpMaxBscan,deltaMin,deltaMax);

[s1,s2,s3,s4] = size(voxels);
layers = findLayers(labels,s1,s2,s3,s4);
RPE1 = cat(1,RPE1,layers{1}+ref-wt);
RPE2 = cat(1,RPE2,layers{2}+ref-wt);

    disp(['end slice ' num2str(i)])

end

%figure;imshow3D(volume_aligned,[],'plot',cat(3,RPE1,RPE2))

end

