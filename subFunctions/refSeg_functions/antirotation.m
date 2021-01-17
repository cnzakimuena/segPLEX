function [ image_out ] = antirotation( image, angle, s1, borderEffect )
if borderEffect
mask = 2*ones(s1(1),s1(2));
mask (1,:) = 1;
mask(s1(1),:) = 1;
mask(:,1) = 1;
mask(:,s1(2))=1;
mask_r = imrotate(mask,angle,'nearest');

    image(mask_r ==0) = 0;
    image(mask_r ==1) = 0;
end

mask_r_ = imrotate(mask_r,-angle,'nearest');
image_out = imrotate(image,-angle,'bilinear');
[row,col] = findRowsCol2Crop(mask_r_);
image_out = image_out(row(1):row(2),col(1):col(2),:);


    
% figure ; imshow(image,[])
%figure; imshow(image_out,[])

end

