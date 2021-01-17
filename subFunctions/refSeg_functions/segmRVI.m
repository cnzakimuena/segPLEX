function [RVI] = segmRVI(volume_aligned,lRVIf,DeltaY)
%SEGMRVI segments the RVI with a 3D graph
%   volume_aligned : 3D OCT volume
%   lRVIf  : first segmentation of the RVI given by the function bscanstore2volume3
%   DeltaY : distance between two adjacent bScan
%   for more details of the 3D graph algorithm,see \readme\Algorithme de segmentation des couches rétiniennes et choroïdiennes.docs) 
%   RVI: results of the segmentation. RVI is a matrix each lines corresponds to the coordinate of the RVI for each Bscan.

%% This function is almost the same as segmRPE. As we already comment segmRPE, we won't comment the following code.
%The main difference between segmRVI and segmRPE is that we don't use the
%inflexion weights but directly the gradient of all pixels. (reminder: the
%inflexion weight of a non inflexion point =0 whereas =gradient(X) in the
%following code)
%There is only one surface to segment whereas in segmRPE there where the
%anterior and posterior 
%%
volume_aligned = single(mat2gray(volume_aligned));
[nbl,nbc,nbz] = size(volume_aligned);
%lRPEmin =  max(min(RPE'));
mask = zeros(nbl,nbc,nbz);
for i =1:nbz
    try
    bscan = zeros(nbl,nbc);
    bscan(sub2ind(size(bscan),min(nbl,max(1,round(lRVIf(i,:)))),1:nbc)) =1;
    mask(:,:,i) = bscan;
    catch
        continue;
    end
end
mask=imdilate(mask,ones(50,1));
[l,~] = find(mask); lmin = min(l); lmax = max(l);
%figure; imshow3D(mask,[],'plot',cat(3,lRVIf))
volRVI =volume_aligned(lmin:lmax,:,:); 
volRVIf = imfilter(volRVI,fspecial('gaussian',[5,5],1),'symmetric');

%figure; imshow3D(volRPEf );
grad =imfilter(volRVIf,[-1;-1;1;1],'symmetric');
gradL1 = grad;%gradL2 = -grad;
gradL1(gradL1<0) = 0; %gradL2(gradL2<0)=0;
gradL1 = gradL1.*mask(lmin:lmax,:,:);
%figure; imshow3D(gradL1,[])
% inflex2 = zeros(size(inflex1));
% inflex2(wt+1,:,:) = 1;
deltaMin=[1,1,1]; deltaMax=[25,25,25];%useless
jumpMaxCol   = 2;
if DeltaY >0.07
    jumpMaxBscan =14;
else
    jumpMaxBscan = 7;
end
%figure; imshow3D(gradL1(:,:,:),[],'plot',cat(3,RPE))

%graph Cut
% if B-scans are larger than design (512x496), graph will be excuted in two
% parts, see '(mod)' below
if size(volume_aligned,2) > 492
    nbBscan = floor(nbz/20); % change dividing number according value of size(volume_aligned,2)
                             % until suitable to available RAM
else
    nbBscan = nbz;
end

RVI = [];
for i =1:floor(nbz/nbBscan)
    
    disp(['begin slice ' num2str(i)])
    
    if length((i-1)*nbBscan+1:min((i+1)*nbBscan,nbz))<nbBscan*2
        inter = (i-1)*nbBscan+1:min((i+1)*nbBscan,nbz);
    else
        inter = (i-1)*nbBscan+1:min((i)*nbBscan,nbz);
    end
    voxels = cat(4,gradL1(:,:,inter)*255);
    voxels (isnan(voxels)) = 0;
    [~, labels] = graphCutPierre2(voxels,jumpMaxCol,jumpMaxBscan,deltaMin,deltaMax);
    [s1,s2,s3,s4] = size(voxels);
    layers = findLayers(labels,s1,s2,s3,s4);
    
    RVI = cat(1,RVI,layers{1}+lmin);
    
    disp(['end slice ' num2str(i)])

end
%figure;imshow3D(voxels,[])
%figure;imshow3D(voxels,[],'plot',cat(3,RVI))
%figure;imshow3D(volume_aligned,[],'plot',cat(3,RVI))
end

