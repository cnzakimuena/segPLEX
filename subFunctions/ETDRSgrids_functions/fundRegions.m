
function regionsETDRS = fundRegions(fundIm, centerXIm, centerYIm, sizeRed)

% figure;imshow(fundIm,[])
% hold on
% plot(centerXIm,centerYIm,'.r')

% fundus image properties ~9x9mm
resFund = 1.95/sizeRed; % fundus resolution um/px
fundDimX = size(fundIm,2)*3*resFund; % actual fundus X dimension in um 
fundDimY = size(fundIm,1)*3*resFund; % actual fundus Y dimension in um 

% *steps*

% (1) foveal location is automatically detected using layers thickness;
% initially searching the lowest thickness in a 1.5x1.5mm circular area, 
% (assuming fundus image is centered at fovea) 
% (2) annular circles are made by assigning 1s to an inner circle, thereby 
% creating an annular mask
% (3) annular quadrants are made by assigning 1s to all areas left or 
% right, top or bottom of the desired quadrants
% (4) angled annular quadrants are made by rotating the mask around the 
% foveal center point 
% (5) a 3D angled annular quadrant mask is obtained by duplicating the 
% initial 2D mask in the axial direction

% ***(OD numenclature-based) ETDRS regions algorithm***

% Create a logical image of a circle with specified diameter, center, and 
% fundus image size

% create the fundus size image ~9x9mm
imageSizeX = round(fundDimX/resFund);
imageSizeY = round(fundDimY/resFund);
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
centerX = centerXIm + round((imageSizeX-round(size(fundIm,2)))/2);
centerY = centerYIm + round((imageSizeY-round(size(fundIm,1)))/2);

% create the circle in the image
diaDim = 3000; % real desired diameter dimension in um 
radiusFac = 1/resFund; % conversion factor
radius = round((diaDim*radiusFac)/2);

% *ETDRS 1-5 —  full inner area mask (3x3mm)*
% circlePixels is a 2D logical array
circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;
%figure; imshow(circlePixels,[])
% hold on
% plot(centerX,centerY,'.r')

% *ETDRS 2-5 —  inner annular area mask (3x3mm)*
diaDim3 = 1000; % real desired diameter dimension in um 
radius3 = round((diaDim3*radiusFac)/2);
circlePixels4 = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius3.^2;
circlePixels5 = circlePixels;
circlePixels5(circlePixels4) = 0;
%figure; imshow(circlePixels5,[])
quadPixels_innerA = circlePixels5;
quadPixels_innerA(:,1:centerX) = 0; 
quadPixels_innerA(centerY:end,:) = 0;
quadPixels_innerA1 = rotateAround(quadPixels_innerA, centerY, centerX, 45); %      ETDRS 2 —   top quadrant 
quadPixels_innerA2 = rotateAround(quadPixels_innerA, centerY, centerX, -45); %     ETDRS 3 —   right quadrant          
quadPixels_innerA3 = rotateAround(quadPixels_innerA, centerY, centerX, -135); %    ETDRS 4 —   bottom quadrant         
quadPixels_innerA4 = rotateAround(quadPixels_innerA, centerY, centerX, 135); %     ETDRS 5 —   left quadrant         
%figure; imshow([quadPixels_innerA1 quadPixels_innerA2 quadPixels_innerA3 quadPixels_innerA4],[])

% *ETDRS 1 —  full center area mask (1x1mm)*
%figure; imshow(circlePixels4,[])

% cropping is based on actual center of the image, not fundus center
imageCenter = round([imageSizeX imageSizeY]/2); %actual image center
imageCenterX = imageCenter(1);
imageCenterY = imageCenter(2);
startX = imageCenterX-round(size(fundIm,2)/2);
endX = startX+size(fundIm,2)-1;
startY = imageCenterY-round(size(fundIm,1)/2);
endY = startY+size(fundIm,1)-1;

% regionsETDRS = zeros([size(fundIm) 9]);
regionsETDRS = zeros([size(startY:endY, 2) size(startX:endX, 2) 5]);
regionsETDRS(:,:,1) = circlePixels4(startY:endY, startX:endX);    
regionsETDRS(:,:,2) = quadPixels_innerA1(startY:endY, startX:endX);       
regionsETDRS(:,:,3) = quadPixels_innerA2(startY:endY, startX:endX);  
regionsETDRS(:,:,4) = quadPixels_innerA3(startY:endY, startX:endX);    
regionsETDRS(:,:,5) = quadPixels_innerA4(startY:endY, startX:endX); 
%figure; imshow(circlePixels4(startY:endY,startX:endX),[])
% hold on
% plot(centerXIm,centerYIm,'.r')
%figure;imshow3D(regionsETDRS,[])

end
