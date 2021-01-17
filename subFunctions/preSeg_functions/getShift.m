
function [sRecord] = getShift(vol, numScan)

sRecord = zeros(numScan,2);

% get downstream 'getRetinaAndBm3' function padding shifting values and
% eliminate their outliers
for i = 1:size(vol,3)
    img = vol(:,:,i);
    %figure; subplot(1,2,1),imshow(img,[]),subplot(1,2,2),imshow(imgf,[])
    imgf = imfilter(img,fspecial('gaussian',[10,10],4),'symmetric'); % gaussian filter
    [~,~, sVals] = findApproxRVIandRPE(imgf); % find a very first approx of the RVI and RPE
    
    % makes column vectors of shifting values on the left and right
    % of the B-scan, i.e., each B-scan in 'findApproxRVIandRPE' above
    % was divided into ten sctions, only the mean RVI locations of the
    % two first sections, and two last sections were stored into
    % 'sVals'
    % 'getRetinaAndBm3' below requires padding which is achieved by
    % duplicating the extreme left and right sections of the B-scan and
    % and concatenating the duplicated sections . The duplicated
    % sections then need to be shifted to account for the desired
    % segmentations' curves
    sRecord(i,:) = [sVals(2)-sVals(1) sVals(3)-sVals(4)];
end
% loop to eliminate outliers
for i = 1:size(sRecord, 1)
    tempVec = [];
    tempVec2 = [];
    if i < (size(sRecord, 1)-10)
        % get middle point of a 10 elements array that follows the current
        % index 'i' within 'sRecord', in which half the numbers are
        % above the median and half are below
        medVal = median(sRecord(i:i+9,1));
        medVal2 = median(sRecord(i:i+9,2));
        % if the current of 10 shifting values from 'sRecord' is less than 5
        % pixels away from the median of a range of 10 elements within
        % which it is found, store it in 'tempVec'
        % the values kept in 'tempVec' are then averaged and assigned to
        % the current 'sRecord' element (of its whole range), thereby having
        % eliminated outliers
        for ii = 1:10
            medDiff = abs(sRecord(i+ii-1,1)-medVal);
            medDiff2 = abs(sRecord(i+ii-1,2)-medVal2);
            if medDiff < 5
                tempVec = [tempVec ; sRecord(i+ii-1,1)];
            end
            if medDiff2 < 5
                tempVec2 = [tempVec2 ; sRecord(i+ii-1,2)];
            end
        end
    elseif i > (size(sRecord, 1)-11)
        medVal = median(sRecord(i-9:i,1));
        medVal2 = median(sRecord(i-9:i,2));
        for ii = 1:10
            medDiff = abs(sRecord(i-ii+1,1)-medVal);
            medDiff2 = abs(sRecord(i-ii+1,2)-medVal2);
            if medDiff < 5
                tempVec = [tempVec ; sRecord(i-ii+1,1)];
            end
            if medDiff2 < 5
                tempVec2 = [tempVec2 ; sRecord(i-ii+1,2)];
            end
        end
    end
    sRecord(i,1) = round(mean(tempVec));
    sRecord(i,2) = round(mean(tempVec2));
end

end