function [ret,bm,yRPEt,yRPEb] = getRetinaAndBm3(im, lShift, rShift, CP)
%GETRETINAANDBM3 : find the position of the retina, bm, rpe (top and
%bottom). 
%   ret : position of the retina after alignment
%   bm : position of the bm before (!) alignment
%   yRPEt : position of the anterior interface RPE after alignment 
%   yRPEb : position of the posterior interface of the RPE after alignment

im4RVI = imopen(im,strel('disk',4,8)); %to segment the RVI we do an opening 
%in order to delete upper noise. 
%figure; imshow(im4RVI,[])

imF = imfilter(im,fspecial('gaussian',[3,3],1),'symmetric'); %thin filter
imFF = imfilter(im,fspecial('gaussian',[6,6],3),'symmetric');%large filter
imFF2 = imFF;
imFF = imclose(imFF,strel('disk',2,8));%closing is to uniformise the RPE  
imFF2 = imclose(imFF2,strel('disk',2,8)); 
% imF(imF>prctile(imF(:),99.5)) =prctile(imF(:),99.5);
% imFF(imFF>prctile(imFF(:),99.95)) =prctile(imFF(:),99.95);
sz = size(imF); %size
padSize = round(size(imF,2)/10); %nb of column to add on the left and on the right 

%in order to rely on Dijkstra algorithm on the border of the image
lpadFac = abs(2.1992*exp(-0.012*lShift)); % fitted correction factor
rpadFac = abs(2.1992*exp(-0.012*rShift));
imF = [imtranslate(imF(:,1:padSize),[0, -lShift*lpadFac]) imF imtranslate(imF(:,(end-padSize+1):end),[0, -rShift*rpadFac])]; 
imFF = [imtranslate(imFF(:,1:padSize),[0, -lShift*lpadFac]) imFF imtranslate(imFF(:,(end-padSize+1):end),[0, -rShift*rpadFac])];
imFF2 = [imtranslate(imFF2(:,1:padSize),[0, -lShift*lpadFac]) imFF2 imtranslate(imFF2(:,(end-padSize+1):end),[0, -rShift*rpadFac])];
im4RVI = [imtranslate(im4RVI(:,1:padSize),[0, -lShift*lpadFac]) im4RVI imtranslate(im4RVI(:,(end-padSize+1):end),[0, -rShift*rpadFac])];
% % imshow(imF,[])

%imFF=imclose(imFF,strel('disk',12,8));
%figure; imshow(im,[])
%figure; imshow(imF,[])
%figure; imshow(imFF,[])
sizFil = 2;
%imStep    = imfilter(imFF, heaviside(-sizFil:sizFil)' - 0.5);  
imStep    = imfilter(imFF, [-1;0;0;1.08]);  %axial gradient (ascendent) [-0.9;0;0;1.2]

imStep2 = imfilter(imFF2, [-1;0;0;1]);

imStepRVI    = imfilter(im4RVI, [-1;0;0;1]); %axial gradient(ascendent)
imStepInv = - imStep2; % axial gradient (descendent)
%figure; subplot(2,1,1) ,imshow(imFF,[]),subplot(2,1,2),imshow(imStepInv,[])
imStepInv(imStepInv < 0) = 0; %delete negative value
imStepInv([1:2*sizFil,end-2*sizFil:end],:) = 0; %delete information in the up and bottom rows
%figure; imshow(imStepInv,[])

% Reduce smoothing to get better precision
imStepPrecise = - imfilter(imF, [-1;1]);%precise gradient
imStepPrecise(imStepPrecise < 0) = 0;%delete negative values
imStepPrecise([1:2*sizFil,end-2*sizFil:end],:) = 0;%delete information in the up and bottom rows
% imStepPreciseBW = imbinarize(imStepPrecise);
% imStepPreciseBW = imclose(imStepPreciseBW,strel('disk',1,8));
% imStepPreciseBW = bwareaopen(imStepPreciseBW, 10);
% imStepPrecise = imStepPreciseBW.*imStepPrecise;

yFirst  = NaN(1,size(imStep,2));%location of the first positive peak =RVI
yRPEt = NaN(1,size(imStep,2));%location of the second positive peak =RPEt
yThird  = NaN(1,size(imStep,2));%location of the first negative peak = RPEb 
%%

% loop below finds filtered image maximum intensities, and selects
% the values of the pixels at indexes nearer the top locations of the
% filtered image, which are representative of the RVI
for k = 1:size(imStep,2) %for each bscan, segment the RVI 
    imStepRVI_i = imStepRVI(:,k);
    [pks,locs] = findpeaks(imStepRVI_i,'SortStr','descend','MinPeakDistance',15);
    ind = pks<max(pks)/4;%delete peaks which intensity is lower than a fourth of the main peak
    pks2 =pks; locs2 = locs;
    pks2(ind)=[]; 
    locs2(ind)=[];
    if numel(locs2) < 2, continue, end
    [locs2, ix] = sort(locs2);%sort the location 
    pks2        = pks2(ix);
    yFirst(k)  = locs2(1);
    
end
yFirst  = eliminateJumps(yFirst);

%figure;imshow(im4RVI,[]),hold on, plot(yFirst)
%%

for k = 1:size(imStep,2)%for each ascan, segment the top interface of the RPE

    imStep_i = imStep(:,k);
    [pks,locs] = findpeaks(imStep_i,'SortStr','descend','MinPeakDistance',15);
    if numel(locs) < 2, continue, end
    ind2Del = locs<yFirst(k)+20; %to delete peaks before the RVI
    pks(ind2Del) = []; locs(ind2Del) =[]; 
    if numel(locs) < 2, continue, end
    [~,indM] = max(pks);
    yRPEt(k) = locs(indM);
    
end

yRPEt = eliminateJumps(yRPEt);%refine the segmentation by eliminating mistakes
%the refinment is done with Dijkstra algorithm ( = shortest path algorithm)

% figure;imshow(imF,[]),hold on, plot(yRPEt), plot(yFirst)

for  k = 1:size(imStep,2)%for each scan, segment the bottom of the RPE
    [pks,locs] = findpeaks(imStepInv(:,k),'MinPeakDistance',7);
    % based on change in axial resolution, assume PLEX RPE thickness 2x 
    % that observed using Spectralis ('yRPEt(k)+21' becomes 'yRPEt(k)+42')
    ind2Del = locs>yRPEt(k)+42|locs < yRPEt(k)+10; 
    locs( ind2Del) = [];
    pks( ind2Del) = [];
    ind2Del = pks<max(pks)/2;
    locs( ind2Del) = [];
    pks( ind2Del) = [];
    [~,ix] = max(pks(1:min(length(pks),2)));%find(locs > ySecond(k),1,'first');
    if isempty(ix), continue, end
    
    yThird(k) = locs(ix);
    %figure; plot( imStep(:,k)),hold on, plot(locs,-pks,'xr')
end

yThird  = eliminateJumps(yThird);%refine the segmentation by eliminating mistakes
%the refinment is done with Dijkstra algorithm ( = shortest path algorithm)
yThird(yThird<yRPEt) =yRPEt(yThird<yRPEt);
%figure;imshow(imFF,[]),hold on, plot(ySecond), plot(yFirst), plot(yThird)
%figure;imshow(imF,[]),hold on,plot(yFirst), plot(ySecond),plot(yThird)
%figure;imshow(imF,[]),hold on,plot(yFirst), plot(yRPEt),plot(yThird)

[yRPEb] = findRPEbottom2(imF,imStepPrecise,yRPEt, yThird);%refine the RPE bottom segmentation

%figure;imshow(imF,[]),hold on, plot(yFirst),plot(yRPEt),plot(yRPEb,'x'),plot(yThird)

%delete padded columns
yFirst(1:padSize) = [];yFirst(end-padSize+1:end)=[];
yRPEt(1:padSize) = [];yRPEt(end-padSize+1:end)=[];
yThird(1:padSize) = [];yThird(end-padSize+1:end)=[];
yRPEb(1:padSize) = [];yRPEb(end-padSize+1:end)=[]; 

%figure;imshow(im,[]),hold on, plot(yFirst),plot(yRPEt),plot(yRPEb,'x'),plot(yThird)

% Compute a convex-hull below the RPE bottom limit to estimate the Bruch's
% membrane. Don't really understand this part but it works ! (Pierre)
convexEnd = sz(2)-20;
DT1 = DelaunayTri([1:convexEnd]',yRPEb(1:convexEnd)');
CH1 = convexHull(DT1);
CHpts1 = flipud([DT1.X(CH1,1) DT1.X(CH1,2)]);
last1 = find(CHpts1(:,1) < circshift(CHpts1(:,1),1),1,'first') - 1;
CHcurve1 = fit(CHpts1(1:last1,1),CHpts1(1:last1,2),'linear');
CHcurve1 = round(CHcurve1(1:convexEnd));

DT2 = DelaunayTri([convexEnd+1:length(yRPEb)]',round(yRPEb(convexEnd+1:end))');
%DT2 = delaunayTriangulation([convexEnd+1:length(yRPEb)]',round(yRPEb(convexEnd+1:end))');

if~isempty(DT2.Triangulation)
CH2 = convexHull(DT2);

CHpts2 = flipud([DT2.X(CH2,1) DT2.X(CH2,2)]);
last2 = find(CHpts2(:,1) < circshift(CHpts2(:,1),1),1,'first') - 1;
CHcurve2 = fit(CHpts2(1:last2,1),CHpts2(1:last2,2),'linear');
CHcurve2 = round(CHcurve2(convexEnd+1:sz(2)));
else
    CHcurve2 = yRPEb(convexEnd+1:end)';
end
CHcurve=[CHcurve1;CHcurve2];
% Set results
ret = yFirst;  % First retina layer
bm  = CHcurve; % Estimation of Bruch's membrane
%figure;imshow(im,[]),hold on, plot(yRPEt),plot(yThird), plot(yFirst)
%figure;imshow(im,[]),hold on, plot(yThird), plot(bm);

%figure; imshow(imF)

end

function [x,y] = traceGraph(aC,bC,aIm,bIm,imind,num,edges,bscan)

bscan = mat2gray(double(bscan));
[m,n] = size(bscan);

startedge=edges(:,1);
endedge=edges(:,2);
startlength=length(find(startedge));
endlength=length(find(endedge));

% DL = - imfilter(bscan,[-1;1],'symmetric');
DL = bscan;
DL(DL<0)=0;
DL=mat2gray(DL);

s=zeros(size(aC));
s(aC~=1 & bC~=num) = 2 - (DL(aIm) + DL(bIm));
s(1:startlength)=1;
s(length(s)-endlength+1:length(s))=1;
% Each element in s correspond to a graph edge
C=sparse(aC,bC,s,num,num);
[~,path,~]=graphshortestpath(C,1,num); %path is the index in the nodes array
[y,x]=ind2sub([m,n],imind(path(2:end-1)-1));

check=[[x;0] [0;x]];
check=(check(:,1)==check(:,2));
check=find(check);
y=y(setdiff(1:length(x),check));
x=x(setdiff(1:length(x),check));

end

function [msk,outWeigth] = makeRPEmask(imOri,inWeigth,top, bottom, retina)

% Creates a mask with foreground pixels in a strip between "top" and
% "bottom" traces. It fills up the columns that do not have data for the
% traces with the nearest column information

sz = size(inWeigth);

msk = zeros(sz);

outWeigth = inWeigth;

% Find the start and end of all the traces together
start = find(~isnan(top) & ~isnan(bottom) & ~isnan(retina),1,'first');
fin   = find(~isnan(top) & ~isnan(bottom) & ~isnan(retina),1,'last');

rpeThickness = nanmedian(bottom - top);

% Computes the absolut top limit of the region (halfway between retina and RPE)
% gapTop = nanmedian(top - retina) / 2;
% topLim = round(min(top, retina + gapTop));
topLim = NaN(size(top));

for k = start:fin
      topLim(k) = round(sum(imOri(top(k):bottom(k),k) .* (top(k):bottom(k))') / sum(imOri(top(k):bottom(k),k)));
end

% Computes the absolut bottom limit of the region (halfway between RPE and bottom of image)
% gapBot = nanmedian(sz(1) - bottom) / 2;
% botLim = round(max(sz(1) - gapBot, bottom));
botLim = round(min(sz(1), bottom + rpeThickness / 2));

for k = start:fin
      msk(topLim(k):botLim(k),k) = true;
%     msk(topLim(k):botLim(k),k) = true;
%     outWeigth([1:top(k),bottom(k):end],k) = 0.1 * outWeigth([1:top(k),bottom(k):end],k);
end

% Fill up gaps at start and end
if start ~= 1
    msk(topLim(start):botLim(start),1  :start) = true; 
    
%     rng = [1:top(start),bottom(start):sz(1)];
%     outWeigth(rng,1:start) = 0.1 * outWeigth(rng,1:start);
end


if fin ~= sz(2)
    msk(topLim(fin):botLim(fin),fin:end)   = true; 
    
%     rng = [1:top(fin),bottom(fin):sz(1)];
%     outWeigth(rng,fin:end) = 0.1 * outWeigth(rng,fin:end);
    
end

outWeigth(~msk) = 0.1 * outWeigth(~msk); 

end

function traceOut2 = eliminateJumps(traceIn)
%ELIMINATEJUMPS eliminates potential jumps in the first rough segmentation
%   TRACEIN is the position of the first segmentation
%   TRACEOUT2 is the position of the segmentation after correction using
%   Dijkstra algorithm

xMax = length(traceIn);
nodeMask = zeros(max(traceIn),length(traceIn));%mask of the nodes in the image
x= 1:length(traceIn);
x(isnan(traceIn))=[];
traceIn(isnan(traceIn)) =[];

nodeMask(sub2ind(size(nodeMask),round(traceIn),x))=1;
[~,traceOut2] = graphSearch2(nodeMask,'graph1');%connect nodes eliminating errors using Dijkstra algorithm
traceOut2 = interp1(traceOut2.X,traceOut2.Y,1:xMax,'linear'); % interpolate the filtered nodes to obtain a smooth line
%figure; imshow(nodeMask)

%Extrapolate the begining and the and using the nearest neighbour
%(only useful in case of NAN at the beginning or at the end) 
start = find(~isnan(traceOut2),1,'first');
fin   = find(~isnan(traceOut2),1,'last');

traceOut2(1:start) = traceOut2(start);
traceOut2(fin+1:end) = traceOut2(fin);

end

function [xOut,traceOut] = findRPEbottom(imOri,invGrad,top, bottom)

sz = size(invGrad);

msk = false(sz);

% Find the start and end of all the traces together
start = find(~isnan(top) & ~isnan(bottom),1,'first');
fin   = find(~isnan(top) & ~isnan(bottom),1,'last');

top(1:start) = top(start);
top(fin:end) = top(fin);
bottom(1:start) = bottom(start);
bottom(fin:end) = bottom(fin);
% Computes RPE msk
%figure; imshow(imOri,[]),hold on, plot(top),plot(bottom)
for k = 1:numel(top)
    if round(top(k))~=round(bottom(k))
        if round(bottom(k))-round(top(k))>10
            m = (round(bottom(k))-round(top(k)))/2;
            msk(round(top(k)+m):round(bottom(k)),k) = true; %Pierre:-1
        else
            msk(round(top(k)):round(bottom(k)),k) = true; %Pierre:-1
        end
    else
        bande =3;
        m = max(1,round(top(k))-bande);
        M = min(sz(1),round(bottom(k))+bande);
        msk(round(top(k)):round(bottom(k)),k) = true;
    end
end


wInt = imOri(msk);
wInt = (wInt - min(wInt)) / (max(wInt) - min(wInt)); 
wMatrix = zeros(size(msk));
wMatrix(msk) = wInt;

[aC,bC,aIm,bIm,imind,edges,num] = ConnectivityMatrix(msk,8);%figure; imshow(msk,[])
%figure; subplot(2,1,1),imshowpair(msk,imOri),hold on, plot(bottom);subplot(2,1,2),imshow(msk.*imOri)
[xRPE,yRPE] = traceGraph(aC,bC,aIm,bIm,imind,num,edges,wMatrix);%????
cpt = 0;
while(isempty(xRPE)&& cpt <= 5)
    msk = imdilate(msk,[1;1;1;1;1;1;]);
    [aC,bC,aIm,bIm,imind,edges,num] = ConnectivityMatrix(msk,8);%figure; imshow(msk,[])
%figure; subplot(2,1,1),imshowpair(msk,imOri);subplot(2,1,2),imshow(msk.*invGrad)
    [xRPE,yRPE] = traceGraph(aC,bC,aIm,bIm,imind,num,edges,wMatrix);%????
    cpt = cpt+1;
end
xRPE = round(xRPE);
yRPE = round(yRPE);
%figure; imshow(imOri,[]);hold on, plot(xRPE,yRPE)
ix = sub2ind(size(imOri),yRPE,xRPE);
wRPE = imOri(ix);

% Refine BM
yRPE = eliminateJumps(yRPE');

% Estimates the bottom edge of the RPE
rpeThickness = nanmedian(bottom - top);

botLim = round(min(sz(1), yRPE + rpeThickness));
msk(:) = false;

for k = 1:numel(yRPE)
      msk(round(yRPE(k)):round(botLim(k)),k) = true;
end
%figure; imshow(msk,[]);hold on, plot(xRPE,yRPE)
wGrad = invGrad(msk);
wGrad = (wGrad - min(wGrad)) / (max(wGrad) - min(wGrad)); 

wMatrix = zeros(size(msk));
wMatrix(msk) = wGrad;

[aC,bC,aIm,bIm,imind,edges,num] = ConnectivityMatrix(msk,8);
[xBot,yBot] = traceGraph(aC,bC,aIm,bIm,imind,num,edges,wMatrix);
while(isempty(xBot)&& cpt <= 5)
    msk = imdilate(msk,[1;1;1;1;1]);
    [aC,bC,aIm,bIm,imind,edges,num] = ConnectivityMatrix(msk,8);%figure; imshow(msk,[])
%figure; subplot(2,1,1),imshowpair(msk,imOri);subplot(2,1,2),imshow(msk.*imOri)
    [xBot,yBot] = traceGraph(aC,bC,aIm,bIm,imind,num,edges,wMatrix);%????
    cpt = cpt+1;
end
% Estimated distance to BM

estDistBM = nanmedian(yBot' - yRPE(xBot));
xOut = 1:numel(yRPE);
traceOut = yRPE + estDistBM;
%figure; imshow(imOri,[]),hold on, plot(traceOut),plot(bottom)





end
function [yBot2] = findRPEbottom2(imOri,invGrad,top, bottom)

sz = size(invGrad);

msk = false(sz);

% Find the start and end of all the traces together
start = find(~isnan(top) & ~isnan(bottom),1,'first');
fin   = find(~isnan(top) & ~isnan(bottom),1,'last');

top(1:start) = top(start);
top(fin:end) = top(fin);
bottom(1:start) = bottom(start);
bottom(fin:end) = bottom(fin);
% Computes RPE msk
%figure; imshow(imOri,[]),hold on, plot(top),plot(bottom)
for k = 1:numel(top)
    if round(top(k))~=round(bottom(k))
        if round(bottom(k))-round(top(k))>10
            m = floor((round(bottom(k))-round(top(k)))/2);
            msk(max(round(top(k)+m),1):min(round(bottom(k)+1),sz(1)),k) = true; %Pierre:-1
        else
            msk(max(round(top(k)),1):min(round(bottom(k)+1),sz(1)),k) = true; %Pierre:-1
        end
    else
        bande =3;
        m = max(1,round(top(k))-bande);
        M = min(sz(1),round(bottom(k))+bande);
        msk(round(top(k)):round(bottom(k)),k) = true;
    end
end



wGrad = invGrad(msk);
wGrad = (wGrad - min(wGrad)) / (max(wGrad) - min(wGrad)); 

wMatrix = zeros(size(msk));
wMatrix(msk) = wGrad;

[aC,bC,aIm,bIm,imind,edges,num] = ConnectivityMatrix(msk,8);
[xBot,yBot] = traceGraph(aC,bC,aIm,bIm,imind,num,edges,wMatrix);
cpt = 0;
while(isempty(xBot)&& cpt <= 5)
    msk = imdilate(msk,[1;1;1;1;1]);
    [aC,bC,aIm,bIm,imind,edges,num] = ConnectivityMatrix(msk,8);%figure; imshow(msk,[])
%figure; subplot(2,1,1),imshowpair(msk,imOri);subplot(2,1,2),imshow(msk.*imOri)
    [xBot,yBot] = traceGraph(aC,bC,aIm,bIm,imind,num,edges,wMatrix);%????
    cpt = cpt+1;
end
yBot2 = interp1(xBot,yBot,1:size(imOri,2),'linear');


end