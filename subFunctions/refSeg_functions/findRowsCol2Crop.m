function [L,C] = findRowsCol2Crop(image)
[nbL,nbC] = size(image);
image2 = [zeros(1,nbC);image;zeros(1,nbC)];
image2 = [zeros(nbL+2,1),image2,zeros(nbL+2,1)];
hL = [1;-1];
imageL = imfilter(image2,hL);
lines = sum(imageL,2);
[~,indm] = min(lines);
[~,indM] = max(lines);
L = [indm,indM-1];

hC = [1,-1];
imageC = imfilter(image2,hC);
col= sum(imageC,1);
[~,indm] = min(col);
[~,indM] = max(col);
C = [indm,indM-1];
%figure; imshow(imageL,[])
%figure; imshow(image,[])

%figure; imshow(image,[])
end