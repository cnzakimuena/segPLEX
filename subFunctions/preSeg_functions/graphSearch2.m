function [PathPts,usedNodes] = graphSearch2(varargin)
%%GRAPHSEARCH2 use Dijkstra algorithm to eliminate errors in a segmentation
%   VARARGIN{1} : mask of the nodes
%   VARARGIN{2} : type of method used 
%   USEDNODES : coordinates of the nodes selected
%   PathPts : I don't know what this is but it's useless (Pierre G)
    nodesMask = varargin{1};
if nargin>1
    whatGraph = varargin{2};
end
switch whatGraph
    case 'graph1' %first method 
%% graph parameter definition        
delColmax=100;% the graph connections will only exist between nodes which distance is less than 'delColmax' columns
delRowmax=50; % the graph connections will only exist between nodes which distance is less than 'delRowmax' rows   
maxJumpCol=2; % the penality will be applied only if the jump is more than 'maxJumpCol' columns     
maxJumpRow=6; % the penality will be applied only if the jump is more than 'maxJumpRow' rows    
alpha_vert = 20;%vertical penality coefficient    
alpha_hori = 10; %horizontal penality coefficient    
[nRows,nCols] = size(nodesMask);

% Index of nodes
[rows, cols] = find(logical(nodesMask));
% Arrays holding the row and col
rows = [0; rows; 0]; % we add nodes before and after the first and the last nodes 
cols = [0;cols;0];% we add nodes before and after the first and the last nodes 
numNodes = numel(find(logical(nodesMask)))+2; %'+2' because we add nodes before and after
connectMatrix = sparse(numNodes,numNodes);%element ij ==> weight of edge between nodes i and nodes j

nodesMask = bwlabel(logical(nodesMask)); %%cannot explain why this is useful
nodesMask(logical(nodesMask)) = nodesMask(logical(nodesMask)) + 1;

nodesIdx = 1:numNodes;
firstCol   = cols(2);
lastCol    = cols(end-1);

% counts = zeros(1,5);
%% Boundary Conditions (calculate weight between the first nodes and the others and
                                        %between the last nodes and the others 
BoundaryColCon = round(delColmax);
LBoundIndx     = nodesIdx(ismember(cols,firstCol:firstCol+BoundaryColCon));
RBoundIndx     = nodesIdx(ismember(cols,max(1,lastCol-BoundaryColCon):lastCol));

%between the first nodes and the next:
connectMatrix(1,LBoundIndx) = alpha_vert *costFunction(cols(LBoundIndx),maxJumpCol,2,'nonLinear');
%between the last nodes and the previous:
connectMatrix(RBoundIndx,numNodes) =  alpha_vert *costFunction((nCols+1 - cols(RBoundIndx)),maxJumpCol,2,'nonLinear');
%% Calculus of the weight for all other edge
for indx = 2:numNodes-1
    if cols(indx) == max(cols)
        break
    end
    connected = nodesIdx(abs(cols(indx)-cols) <= delColmax  &...
                               cols >  cols(indx) &...
                         abs(rows(indx)-rows) <= delRowmax  &...
                                         cols <= nCols           );
    
    if isempty(connected) && ~any(any(connectMatrix(1:indx-1,indx+1:end)))
        nextcol   = min(cols((cols>cols(indx) & abs(rows(indx)-rows)<=delRowmax)));
        
        if isempty(nextcol), break, end % JM if there are no more colums it stops the loop
        
        connected = nodesIdx(cols==nextcol & abs(rows(indx)-rows)<=delRowmax);
    end
    
    dely=abs(rows(indx)-rows(connected));
    delx=abs(cols(indx)-cols(connected));
    
    VertJumpPenalty   =costFunction(dely,maxJumpRow,2,'nonLinear');
    HorizJumpPenalty  =costFunction(delx,maxJumpCol,2,'nonLinear');
    Weight = alpha_vert*VertJumpPenalty + alpha_hori*HorizJumpPenalty ;
    Weight(Weight<= 3 * eps(Weight)) = 0;
    connectMatrix(indx,connected) = Weight;
    
end

[dist,path,~] = graphshortestpath(connectMatrix,1,numNodes,'directed',true);

% If there is no path it tries to solve in parts.
%% Pierre : 'I don't understand this part of the code' (21/11/2017)
if isempty(path) || any(isinf(dist))%|| length(path)<20
    paths = getPossiblePath(connectMatrix);
    
    if isempty(paths)
        PathPts   = [];
        usedNodes = [];        return
    end
    
    PathPts = [];
    
    %Contains the domains for each stretch candidate
    domains = zeros(numel(paths),size(nodesMask,2)); 
    
    for k = 1:numel(paths)
        
        stretch.x = cols(paths(k).ix);
        stretch.y = rows(paths(k).ix);
        
        stretch.weight     = nodesMask(sub2ind(size(nodesMask),stretch.y,stretch.x));
        stretch.sumWeight  = sum(full(stretch.weight));
        stretch.meanWeight = mean(full(stretch.weight));
        stretch.meanHeight = mean(full(stretch.y));
        stretch.length     = numel(full(stretch.y));
        stretch.keep       = 0;
        
%         % Skip unimportant stretches
%         if stretch.meanWeight < 0.5 || stretch.length < 5
%             stretch.x      = 0;
%             stretch.y      = 0;
%             stretch.weight = 0;
%         else
           %Set ones for each columns where the stretch exists
           domains(k,min(stretch.x):max(stretch.x)) = 1;         
%         end
        
        PathPts = [PathPts stretch];
         
    end
    
    if size(domains,1) > 1
      count = sum(domains);
    else
      count = domains > 0;  
    end
    
    ids   = bi2de(domains')'; % Builds a unique number for each combination of segments
    ix    = setdiff(unique(ids),[0, 2.^(0:size(domains,1)-1)]); % Gets a list of the unique numbers except those representing only one segement
    
    keep = ones(1,size(domains,1));
    
    for k = 1:numel(ix)
          
          objIx = find(de2bi(ix(k),numel(paths))); %List segments that overlap here
          
          objIx = setdiff(objIx,find(keep==0));
          
          sumWeights  = [PathPts(objIx).sumWeight];
          meanWeights = [PathPts(objIx).meanWeight];
          heights     = [PathPts(objIx).meanHeight]; % disfavors segments at the bottom of the scan (originated on noise)
          
          sumWeightRates  = parameters.segmentSelectionSumWeigth  *  sumWeights  / sum(sumWeights);
          meanWeightRates = parameters.segmentSelectionMeanWeigth *  meanWeights / sum(meanWeights);
          heightRates     = parameters.segmentSelectionHeights    ./ heights     / sum(1./heights);
          
          fullRate = sumWeightRates + meanWeightRates + heightRates;
          
%           [~,ixWeight] = sort([PathPts(objIx).meanWeight],'ascend');
%           [~,iHeight]  = sort([PathPts(objIx).meanHeight],'descend');
%           [~,iLength]  = sort([PathPts(objIx).length],'ascend');
          
          % Combine the three ordering criteria
%           [~,ixWin] = max(mean([ixWeight;iHeight;iLength]));
          [~,ixWin] = max(fullRate);
          
          keep(setdiff(objIx,objIx(ixWin))) = 0; %Discard those not choosen 
          
%        end
    end
     
    % Keep segments that do not overlap and remove unsignificant ones
    for k = 1:numel(PathPts)
        if PathPts(k).meanWeight < parameters.minMeanPathWeight || PathPts(k).length < parameters.minSumPathWeigth
            PathPts(k).keep = 0;
        elseif all(count(logical(domains(k,:))) == 1)
            PathPts(k).keep = 1;
        else
            PathPts(k).keep = keep(k);
        end
    end
    
    usedNodes = [];
    
else
  
    usedNodes.X = cols(path(2:end-1));
    usedNodes.Y = rows(path(2:end-1));
    
    PathPts.x      = usedNodes.X;
    PathPts.y      = usedNodes.Y;
    
end


%% GRAPH 3 ____________________________________________________________________
%______________________________________________________________________________
    case 'graph3'

%% THIS PART OF THE CODE IS THE SAME AS THE ONE IN CASE 'GRAPH1' BUT WITH DIFFERENT PARAMETERS. THAT'S WHY IT IS NOT DOCUMENTED AGAIN
parameters = loadParameters;
delColmax=100; %for the connections
delRowmax=100; %for the connections
maxJumpCol=0;% %for the weights
maxJumpRow=0;% % for the weights
alpha = 2;
wM = 20000;
alpha_vert = 1;
alpha_hori = 1;

nbPad = round(size(nodesMask,2)/10);
nodesMask = padarray(nodesMask,[0,nbPad],'symmetric');
[nRows,nCols] = size(nodesMask);

% Index of nodes
[rows, cols] = find(logical(nodesMask));
% Arrays holding the row and col
rows = [0; rows; 0];
cols = [0;cols;0];
numNodes = numel(find(logical(nodesMask)))+2;
connectMatrix = sparse(numNodes,numNodes);

nodesMask = bwlabel(logical(nodesMask));
nodesMask(logical(nodesMask)) = nodesMask(logical(nodesMask)) + 1;

nodesIdx = 1:numNodes;
firstCol   = cols(2);
lastCol    = cols(end-1);

% counts = zeros(1,5);
%% Boundary Conditions
BoundaryColCon = round(delColmax);
LBoundIndx     = nodesIdx(ismember(cols,firstCol:firstCol+BoundaryColCon));
RBoundIndx     = nodesIdx(ismember(cols,max(1,lastCol-BoundaryColCon):lastCol));

connectMatrix(1,LBoundIndx) = 1 + (cols(LBoundIndx).^2+...
    wM  * (heaviside(cols(LBoundIndx)-maxJumpCol)           .*...
    abs((cols(LBoundIndx)-maxJumpCol)))                     .*...
    (sigmf(cols(LBoundIndx),[alpha,maxJumpCol])-0.5)              );
%connectMatrix(1,LBoundIndx) = alpha_hori *costFunction(cols(LBoundIndx),maxJumpCol);
connectMatrix(RBoundIndx,numNodes) = 1 + ((nCols+1 - cols(RBoundIndx)).^2+...
    wM  * (heaviside((nCols+1 - cols(RBoundIndx)) - maxJumpCol)                .*...
    abs(((nCols+1 - cols(RBoundIndx)) - maxJumpCol)))                          .*...
    (sigmf((nCols+1 - cols(RBoundIndx)), [alpha,maxJumpCol])-0.5)                    );
%connectMatrix(RBoundIndx,numNodes) =  alpha_hori *costFunction((nCols+1 - cols(RBoundIndx)),maxJumpCol);
%% All Other Connections
for indx = 2:numNodes-1
    
    
    if cols(indx) == max(cols)
        break
    end
    
    connected = nodesIdx(abs(cols(indx)-cols) <= delColmax  &...
                               cols >  cols(indx) &...
                         abs(rows(indx)-rows) <= delRowmax  &...
                                         cols <= nCols           );
    
    if isempty(connected) && ~any(any(connectMatrix(1:indx-1,indx+1:end)))
        nextcol   = min(cols((cols>cols(indx) & abs(rows(indx)-rows)<=delRowmax)));
        
        if isempty(nextcol), break, end % JM if there are no more colums it stops the loop
        
        connected = nodesIdx(cols==nextcol & abs(rows(indx)-rows)<=delRowmax);
    end
    
    dely=abs(rows(indx)-rows(connected));
    delx=abs(cols(indx)-cols(connected));
    
%     if on4
%         ConnectionAffinity = linePenalty(indx,rows,cols,connected,edgeness);
%     else
%         ConnectionAffinity = Inf;
%     end
%     

%     maxBadPix = 0;
%     nb_badW = mask_badWhite(sub2ind(size(mask_badWhite),rows(connected),cols(connected)));
%     nb_badB = mask_badBlack(sub2ind(size(mask_badWhite),rows(connected),cols(connected)));

   
    %VertJumpPenalty   = (heaviside(dely-maxJumpRow) .* abs((dely-maxJumpRow))) .* (sigmf(dely,[1,maxJumpRow]-0.5));
    VertJumpPenalty   =costFunction(dely,maxJumpRow,2,'nonLinear');
    %HorizJumpPenalty  = (heaviside(delx-maxJumpCol) .* abs((delx-maxJumpCol))) .* (sigmf(delx,[1,maxJumpCol]-0.5));  %was wm/2
    HorizJumpPenalty  =costFunction(delx,maxJumpCol,2,'nonLinear'); 
    %BadWhitePenalty   = wM/2 * (heaviside(nb_badW-maxBadPix).* nb_badW.*sigmf(nb_badW,[alpha,maxBadPix]));
   % BadBlackPenalty   = wM/2 * (heaviside(nb_badB-maxBadPix).* nb_badB.*sigmf(nb_badB,[alpha,maxBadPix]));
    %MEAN_weight =mean([VertJumpPenalty,HorizJumpPenalty ,BadWhitePenalty,BadBlackPenalty]);
    %STD_weight =std([VertJumpPenalty,HorizJumpPenalty ,BadWhitePenalty,BadBlackPenalty]);
    Weight = alpha_vert*VertJumpPenalty + alpha_hori*HorizJumpPenalty ;
%     Weight =    500*VertJumpPenalty         +... % Control of Vertical Jump Penalty Magnitude
%                 BadWhitePenalty +...% Control on the number of White pixel above the indx pixel
%                 BadBlackPenalty ;   % Control ont the number of Black pixel under the indx pixel
%     %Euclid                  +... % Euclidian Distance Squared
    %
        %on1 * HorizJumpPenalty  +... %Control of Horizontal Jump Penalty Magnitude
    Weight(Weight<= 3 * eps(Weight)) = 0;
        connectMatrix(indx,connected) = Weight;
        
end

%Test JM
% connectMatrix = connectMatrix';
%connectMatrix = tril(connectMatrix + connectMatrix');

% connectMatrix(connectMatrix<= 3 * eps(connectMatrix)) = 0;



% [dist,path,~] = graphshortestpath(connectMatrix,1,numNodes);
% test JM
[dist,path,~] = graphshortestpath(connectMatrix,1,numNodes,'directed',true);

% If there is no path it tries to solve in parts.
if isempty(path) || any(isinf(dist))%|| length(path)<20
    paths = getPossiblePath(connectMatrix);
    
    if isempty(paths)
        PathPts   = [];
        usedNodes = [];        return
    end
    
    PathPts = [];
    
    %Contains the domains for each stretch candidate
    domains = zeros(numel(paths),size(nodesMask,2)); 
    
    for k = 1:numel(paths)
        
        stretch.x = cols(paths(k).ix);
        stretch.y = rows(paths(k).ix);
        
        stretch.weight     = nodesMask(sub2ind(size(nodesMask),stretch.y,stretch.x));
        stretch.sumWeight  = sum(full(stretch.weight));
        stretch.meanWeight = mean(full(stretch.weight));
        stretch.meanHeight = mean(full(stretch.y));
        stretch.length     = numel(full(stretch.y));
        stretch.keep       = 0;
        
%         % Skip unimportant stretches
%         if stretch.meanWeight < 0.5 || stretch.length < 5
%             stretch.x      = 0;
%             stretch.y      = 0;
%             stretch.weight = 0;
%         else
           %Set ones for each columns where the stretch exists
           domains(k,min(stretch.x):max(stretch.x)) = 1;         
%         end
        
        PathPts = [PathPts stretch];
         
    end
    
    if size(domains,1) > 1
      count = sum(domains);
    else
      count = domains > 0;  
    end
    
    ids   = bi2de(domains')'; % Builds a unique number for each combination of segments
    ix    = setdiff(unique(ids),[0, 2.^(0:size(domains,1)-1)]); % Gets a list of the unique numbers except those representing only one segement
    
    keep = ones(1,size(domains,1));
    
    for k = 1:numel(ix)
          
          objIx = find(de2bi(ix(k),numel(paths))); %List segments that overlap here
          
          objIx = setdiff(objIx,find(keep==0));
          
          sumWeights  = [PathPts(objIx).sumWeight];
          meanWeights = [PathPts(objIx).meanWeight];
          heights     = [PathPts(objIx).meanHeight]; % disfavors segments at the bottom of the scan (originated on noise)
          
          sumWeightRates  = parameters.segmentSelectionSumWeigth  *  sumWeights  / sum(sumWeights);
          meanWeightRates = parameters.segmentSelectionMeanWeigth *  meanWeights / sum(meanWeights);
          heightRates     = parameters.segmentSelectionHeights    ./ heights     / sum(1./heights);
          
          fullRate = sumWeightRates + meanWeightRates + heightRates;
          
%           [~,ixWeight] = sort([PathPts(objIx).meanWeight],'ascend');
%           [~,iHeight]  = sort([PathPts(objIx).meanHeight],'descend');
%           [~,iLength]  = sort([PathPts(objIx).length],'ascend');
          
          % Combine the three ordering criteria
%           [~,ixWin] = max(mean([ixWeight;iHeight;iLength]));
          [~,ixWin] = max(fullRate);
          
          keep(setdiff(objIx,objIx(ixWin))) = 0; %Discard those not choosen 
          
%        end
    end
     
    % Keep segments that do not overlap and remove unsignificant ones
    for k = 1:numel(PathPts)
        if PathPts(k).meanWeight < parameters.minMeanPathWeight || PathPts(k).length < parameters.minSumPathWeigth
            PathPts(k).keep = 0;
        elseif all(count(logical(domains(k,:))) == 1)
            PathPts(k).keep = 1;
        else
            PathPts(k).keep = keep(k);
        end
    end
    
    usedNodes = [];
    
else
  
    usedNodes.X = cols(path(2:end-1));
    usedNodes.Y = rows(path(2:end-1));
    [x_CSI_graph3,y_CSI_graph3] = deletePad(usedNodes.X,usedNodes.Y,nbPad,nCols);
    PathPts.x      = x_CSI_graph3;
    PathPts.y      = y_CSI_graph3;
%     PathPts.pixWeight = mask_badWhite(sub2ind(size(mask_badWhite),rows(path(2:end-1)),cols(path(2:end-1))))+...
%                      mask_badBlack(sub2ind(size(mask_badWhite),rows(path(2:end-1)),cols(path(2:end-1))));
    %PathPts.keep   = 1;
    
%     vals=fit([0;usedNodes.X;nCols+1],[rows(path(2));usedNodes.Y;rows(path(end-1))],'linear');
%     vals=vals(1:nCols);
%     vals=round(smooth(vals,50,'rloess'))'; %0.2
%     PathPts=sub2ind([nRows nCols],vals,1:length(vals));
    
end



end
