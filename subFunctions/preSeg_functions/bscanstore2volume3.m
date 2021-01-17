
function [vol_f,lRVIf,lBM,lRPEtf,lRPEbf,err,trans_out] = bscanstore2volume3(bscanstore)
%BSCANSTORE2VOLUME3 align all bscans to create a volume. To do that, it
%segments the Bruch Membrane(BM) and shifts all A-scan so that the BM is
%straight (at the 165th line)
%   Input : -BSCANSTORE is a cell containing all bscan (after cropping !)
%   Output :-VOL_F : is a 3D matrix = 3D OCT volume after alignment (to 
%            observe it, use the imshow3D function) 
%           -LRVIF,LRPETF, LRPEbf, are respectively 2D matrix storing the
%           position of the RVI,RPEtop,RPEbottom after alignment
%           -lBM is a 2D matrix containing the position of the BM before(!)
%           alignment for each bscan.
%   Method :all the layers are found by searching the maximim of Ascan
%   gradient. The lines are drawn with the shortest path algorithm
%   (=Dijkstra algorithm).
%   RVI = Vetinal vitreous Interface
%   RPE = Retinal Pigment Epithelium
%   BM = Bruch Membrane
nbScan = numel(bscanstore); % nb of bscan in the volume
first = 1;
% while isempty(bscanstore{first})
%     first=first+1; %find the first non empty bscan (when error occurs, it may happen that some bscan are empty (very rare).)
% end
[nbl,nbc] = size(bscanstore{first});%find the nb of lines and colums of all bscans
volume = single(zeros(nbl,nbc,nbScan));
for i = 1:nbScan
    if~isempty(bscanstore{i})
        volume(:,:,i) = single(mat2gray(bscanstore{i})); %normalisation of the intensity between 0 and 1
    end
end

% figure; imshow3D(volume,[]);

%%

lBM = zeros(nbScan,nbc);
lRVI = zeros(nbScan,nbc);
lRPEt = zeros(nbScan,nbc);
lRPEb = zeros(nbScan,nbc);
err=[];

% get downstream 'getRetinaAndBm3' function padding shifting values and 
% eliminate their outliers
sRec = getShift(volume, nbScan);

% for each bscan, segment the RVI, RPEtop,RPEbottom,and BM.
for i = 1:size(volume,3) 
    
    disp(i)
    try
        img = volume(:,:,i);
        %figure; subplot(1,2,1),imshow(img,[]),subplot(1,2,2),imshow(imgf,[])
        imgf = imfilter(img,fspecial('gaussian',[10,10],4),'symmetric'); % gaussian filter
        [rvi_a,rpe_a, ~] = findApproxRVIandRPE(imgf); % find a very first approx of the RVI and RPE
        rvi_a = max(rvi_a-150,1);% to be sure that min_l is really the minimum axial position
        rpe_a = min(rpe_a+150,size(img,1));%to be sure that max_l is really the max position
        img2 = mat2gray(img(rvi_a:rpe_a,:));%crop the bscan around the RVI and RPE in order to reduce complexity of the analysis
        
        lPad = sRec(i,1);
        rPad = sRec(i,2);

        %figure; imshow(img2,[])
        
        [ret,bmConvex,RPEt,RPEb] = getRetinaAndBm3(img2, lPad, rPad, i);%find position of retina, bm, RPEtop and RPE botom
        
        %the segmentation lines below are for the cropped image, the
        %maximum axial position is added back onto them to place them back
        %into the proper location in reference to the original bscan
        lBM(i,:) = bmConvex+rvi_a;
        lRVI(i,:) = ret+rvi_a;
        lRPEt(i,:) = RPEt+rvi_a;
        lRPEb(i,:) = RPEb+rvi_a;
    catch
        lBM(i,:) = nan;
        lRVI(i,:) = nan;
        lRPEt(i,:) = nan;
        lRPEb(i,:) = nan;
        err = [err,i];
        continue;
    end
%     figure;imshow(volume(:,:,i),[]),hold on,plot(lBM(i,:),'g'),plot(lRVI(i,:),'r'),plot(lRPEt(i,:),'b'),plot(lRPEb(i,:),'y')
%     imshow(volume(:,:,i),[])

end
% figure;imshow3D(volume,[],'plot',cat(3,lBM,lRVI,lRPEt,lRPEb),'LineWidth',2)

% [lBM, rawDiff] = clnupSegm(lBM);

%% align each bscan so that the BM is straight 
ref = 510;%reference for the BM
lRVIf = lRVI+(ref-lBM);
lRPEtf = lRPEt+(ref-lBM);
lRPEbf = lRPEb+(ref-lBM);

[vol_f, trans_out] = flattenVol(volume,lBM,ref,err);
%figure; imshow3D(vol_f,[])
%figure; imshow3D(permute(vol_f,[1 2 3]),[])

end
