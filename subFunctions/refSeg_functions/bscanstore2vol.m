
function vol = bscanstore2vol(store)

first = 1;
nbScan = numel(store); % nb of bscan in the volume
[nbl,nbc] = size(store{first});%find the nb of lines and colums of all bscans
vol = single(zeros(nbl,nbc,nbScan));
for i = 1:nbScan
    if~isempty(store{i})
        vol(:,:,i) = single(mat2gray(store{i})); %normalisation of the intensity between 0 and 1
    end
end

end 
