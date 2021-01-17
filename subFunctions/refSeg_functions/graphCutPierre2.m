function [cut, labels] =graphCutPierre2(voxels,deltaX,deltaZ,deltaMin, deltaMax)
%GRAPHCUTPIERRE2. This function is the 3D Graph segmentation.
%[if you want to rename this function to delete the name of its creator,
%feel free. :)]
%   If you want to deeply understand this function, read the article 'Li et al.
%   Optimal Surface Segmentation in Volumetric Images-A Graph-Theoretic Approach.'
%   you will find it in the \readme subfolder. The most important part of the article is the 3rd one : '' Graph Construction''
%   The following documentation will be based on this article.

%   VOXELS is a 4D matrices. Dimension 1,2,3 are the classic 3D space. You
%   can combine several 3D space if  you want to segment several surfaces.
%   The elements of the VOXELS matrix are the pixel weights of the
%   3D OCT volume . (reminder, the basis of the algorithm is to find the
%   surface whose weight is maximal. The surface weight is defined by the
%   sum of its pixel weights).
%   DELTAX is the maximum jump a surface can do between 2 columns. (it is
%       the first smoothness parameter)
%   DELTAY is the maximum jump a surface can do between 2 Bscans (it is the
%       second smoothness parameter)
%   DELTAMIN is a vector whose elements refer to the minimum pixel distance between the surface.
%       exemple1: DELTAMIN = [2]. The 1st and 2nd surface are separated by 2
%       pixels at least.
%       example2: DELTAMIN = [2,5]. The 1st and 2nd surface are separated by 2
%       pixels at least. The 2nd and 3rd surface are separated by 5
%       pixels at least.
%   DELTAMAX is a vector whose elements refer to the maximum pixel distance between the surfaces.
        %example: DELTAMAX = [30,35]. The 1st and 2nd surface are separated
        %by a maximum of 30 pixels. The 2nd and 3rd surface are separated by a maximum of 35
        %pixels.
%   CUT : weight of the 4D cut surface (useless for our algorithm)
%   LABELS : labels for each pixels. 0 if doesn't belong to any surface. 1
%   if belongs to the first surface, 2 to the second etc. (warning : Label
%   is a vector. To transform it into a coordinate map, see findLayers function.
 
nbSurfaces = size(voxels,4);


infin = inf;
%% etape1 : Translation operation to avoid minimum closed set (see the end of 4.1 in Li et al.)

for i =1:nbSurfaces %for each surface do the translation operation (see the end of 4.1 in Li et al.)
    voxels_i =voxels(:,:,:,i);
    if i ==1
        [nbl,nbc,nbz] = size(voxels_i);
        %translation operation to avoid the minimum closed set is empty.
        %(see the end of 4.1 in Li et al.)
        trans = sum(sum(voxels_i(nbl,:,:)))/(nbc*nbz)+1; voxels_i=voxels_i-trans;
        %transformation (soustraction de la ligne du bas)
        voxels2 = [voxels_i(2:nbl,:,:);zeros(1,nbc,nbz,'single')];voxels_i = voxels_i-voxels2;
        voxels(:,:,:,i) =  voxels_i;
    elseif i==nbSurfaces
        %voxels_i=voxels_i(deltaMin(i-1)+1:nbl,:,:);
        [nbl,nbc,nbz] = size(voxels_i);
        %faire une translation pour éviter que le plus petit ensemble soit vide
        trans = sum(sum(voxels_i(nbl,:,:)))/(nbc*nbz)+1;voxels_i=voxels_i-trans;
        %transformation (soustraction de la ligne du bas)
        voxels2 = [voxels_i(2:nbl,:,:);zeros(1,nbc,nbz,'single')]; voxels_i = voxels_i-voxels2;
        voxels(:,:,:,i) =  voxels_i;
    else
        %voxels_i=voxels_i(deltaMin(i-1)+1:nbl-deltaMin(i),:,:);
        [nbl,nbc,nbz] = size(voxels_i);
        %faire une translation pour éviter que le plus petit ensemble soit vide
        trans = sum(sum(voxels_i(nbl,:,:)))/(nbc*nbz)+1; voxels_i=voxels_i-trans;
        %transformation (soustraction de la ligne du bas)
        voxels2 = [voxels_i(2:nbl,:,:);zeros(1,nbc,nbz,'single')];voxels_i = voxels_i-voxels2;
        voxels(:,:,:,i) =  voxels_i;
    end
    
end
clear voxels_i voxels2; %for memory space

%% etape2  : add the S-T nodes. 
%see 4.2 in Li et al : '' The source S is connected to each node V-
%by a directed arc of cost -w(v); every node v in V+ is
%connected to the sink T by a directed arc of cost w(v).

[nbl,nbc,nbz,nbSurfaces] = size(voxels);
terminalWeights = zeros(nbl*nbc*nbz*nbSurfaces,2,'single');
ind_positiveI = find(voxels>=0);
ind_negativeI = find(voxels<0);
terminalWeights(ind_positiveI,1) = voxels(ind_positiveI); %Sink T
terminalWeights(ind_negativeI,2) = -voxels(ind_negativeI);%Source S
clear ind_positiveI ind_negativeI %for memory space

%% etape3 : edge on the pixel below = Intracolumn arcs Ea
%see 3.1 in Li et al. And more precisely eq(1) and (2).

imageInd = reshape (1:nbl*nbc*nbz*nbSurfaces,[nbl,nbc,nbz,nbSurfaces]);
link1 = cat(5,imageInd(1:nbl-1,:,:,:),imageInd(2:nbl,:,:,:));
[s1,s2,s3,s4,~]=size(link1);
edgeWeights = [reshape(link1,[s1*s2*s3*s4,2]) ,infin* ones(s1*s2*s3*s4,1),zeros(s1*s2*s3*s4,1,'single')];
clear link1
%% etape4 : edge on adjacent columns = Intercolumn arcs Er
%see 3.1 in Li et al. And more precisely eq(3).

%<V(x,y,z),V(x+1,y,max(0,z-deltaX))>
link2_ul = imageInd(1:nbl-deltaX,1:nbc-1,:,:);
link2_dr = imageInd(deltaX+1:nbl,1+1:nbc,:,:);
link21   = cat(5,link2_ul,link2_dr);
[s1,s2,s3,s4,~]=size(link21);
link21   = [reshape(link21,[s1*s2*s3*s4,2,1,1]) ,infin* ones(s1*s2*s3*s4,1),zeros(s1*s2*s3*s4,1,'single')];
edgeWeights = cat(1,edgeWeights,link21);
clear link2_ul link2_dr link21

%<V(x,y,z),V(x-1,y,max(0,z-deltaX))>
link2_ur = imageInd(1:nbl-deltaX,1+1:nbc,:,:);
link2_dl = imageInd(deltaX+1:nbl,1:nbc-1,:,:);
link22   = cat(5,link2_ur,link2_dl);
[s1,s2,s3,s4,~]=size(link22);
link22   = [reshape(link22,[s1*s2*s3*s4,2,1]) ,infin* ones(s1*s2*s3*s4,1),zeros(s1*s2*s3*s4,1,'single')];
edgeWeights = cat(1,edgeWeights,link22);
clear link2_ur link2_ur link22 link22

%<V(x,y,z),V(x,y+1,max(0,z-deltaY))>
link2_ua = imageInd(1:nbl-deltaZ,:,1:nbz-1,:);
link2_dp = imageInd(deltaZ+1:nbl,:,1+1:nbz,:);
link23  = cat(5,link2_ua ,link2_dp);
[s1,s2,s3,s4,~]=size(link23);
link23  = [reshape(link23,[s1*s2*s3*s4,2,1]) ,infin* ones(s1*s2*s3*s4,1),zeros(s1*s2*s3*s4,1,'single')];
edgeWeights = cat(1,edgeWeights,link23);
clear link2_ua link2_dp link23 link23

%<V(x,y,z),V(x,y-11,max(0,z-deltaY))>
link2_up = imageInd(1:nbl-deltaZ,:,1+1:nbz,:);
link2_da = imageInd(deltaZ+1:nbl,:,1:nbz-1,:);
link24  = cat(5,link2_up ,link2_da);
[s1,s2,s3,s4,~]=size(link24);
link24  = [reshape(link24,[s1*s2*s3*s4,2,1]) ,infin* ones(s1*s2*s3*s4,1),zeros(s1*s2*s3*s4,1,'single')];
edgeWeights = cat(1,edgeWeights,link24);
clear link2_up link2_da link24

%% étape5 intersurface arc (Es)
%see 3.2 in Li et al. And more precisely eq(4).

for i =1:nbSurfaces-1
    %<V1(x,y,z), V2(x,y,z-du)>
    link1 = cat(5,imageInd(1:nbl-deltaMax(i),:,:,i),imageInd(deltaMax(i)+1:nbl,:,:,i+1));
    [s1,s2,s3,s4,~]=size(link1);
    linkGraph = [reshape(link1,[s1*s2*s3*s4,2]) ,infin* ones(s1*s2*s3*s4,1),zeros(s1*s2*s3*s4,1,'single')];
    edgeWeights = cat(1,edgeWeights,linkGraph);
    %<V2(x,y,z), V1(x,y,z+du)>
    link1 = cat(5,imageInd(1:nbl-deltaMin(i),:,:,i),imageInd(deltaMin(i)+1:nbl,:,:,i+1));
    [s1,s2,s3,s4,~]=size(link1);
    linkGraph = [reshape(link1,[s1*s2*s3*s4,2]) ,zeros(s1*s2*s3*s4,1,'single'),infin* ones(s1*s2*s3*s4,1)];
    edgeWeights = cat(1,edgeWeights,linkGraph);
    %link base set : <V1(0,0,dl), V2(0,0,0)>
    linkGraph = [imageInd(nbl,1,1,i),imageInd(nbl,1,1,i+1),infin,0];
    edgeWeights = cat(1,edgeWeights,linkGraph);
end

%In the previous part, we build the graph. Now we give this graph to an
%algorithm that will find the surface. (with max flow, min cut algorithm).
[cut, labels] = graphCutMex(double(terminalWeights),double(edgeWeights));

end