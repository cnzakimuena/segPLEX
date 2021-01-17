function [min_l,max_l, val] =findApproxRVIandRPE(img)
%FINDAPPROXRVIANDROE : find the approximate max and min position of the RVI and RPE.
%   IMG is the bscan image
%   min_l : approximate axial min position
%   max_l : approximate axial max position
%
% Method : it averages the intensities in each line of the image and
% find the position of the intensity jump.

    sizFil =5;  %size of the gradient mask (quite rough)
    img2 = img; %to separate input from modified data.
    grad = imfilter(img2,heaviside(-sizFil:sizFil)' - 0.5,'symmetric'); %calculate the gradient
    nb_div = 10; %the image will be divided lateraly in nb_div parts. Each part will be analysed separately to find the position of the intensity jump
    [nl,nc]= size(img2);
    col1 = 1; %counter
    
    i = 1;
    val = [];
    
    nbcol = round(nc/nb_div); %nbcol of one division
    min_l = nl;
    max_l = 1;
    
    %for each part of the image, find roughly the position of the RVI
    while(col1<nc) 
        
        % obtain a section of the B-Scan
        grad_i= grad(:,col1:min(col1+nbcol,nc)); 
        
        % delete negative gradient value since the RVI intensities range 
        % from low to high gradient intensity
        grad_i(grad_i<0)=0;  
        
        % vector array of average intensity for each 1536 lines of the 
        % b-scan; 'moy' will have low values except at high gradient 
        % locations
        moy = mean(grad_i,2); 
        
        % intensities and locations (along 'moy' A-scan equivalent array) 
        % of the local maximas; maximas must be at least 11 pixels from
        % each other along the array
        [pks,locs] = findpeaks(double(moy),'SortStr','descend','MinPeakDistance',11);
        
        % obtain only locations of the 1st and 2nd highest intensities, 
        % which should correspond to RVI and RPE (backscatter intensities 
        % along the A-scan direction decrease with tissue depth)
        [~,indm] = min(locs(1:2)); % returns array index (range 1-2)
        [~,indM] = max(locs(1:2)); % returns array index (range 1-2)
        min_l_i = locs(indm);
        max_l_i = locs(indM);
        
        % store the location of 1st and 2nd, 2nd last and last B-scan 
        % sections' RVI locations into a vector array
        colArr = 1:nb_div;
        if i == colArr(1) || i == colArr(2) || i == colArr(end-1) || i == colArr(end) 
            val = [val min_l_i];     
        end       
         
        if min_l_i < min_l % if the jump is before the previous min, then store it
            min_l = min_l_i;
        end
        if max_l_i > max_l % if the jump is after the previous max, then store it
            max_l = max_l_i;
        end
         
        i = i + 1; 
        col1 = col1 + nbcol+1;
    end
    %figure; plot(moy),hold on, plot(locs,moy(locs),'xr')
    %figure; imshow(img2),hold on, plot([1,nc],[min_l,min_l],'g'),hold on, plot([1,nc],[max_l,max_l],'r')
    %figure; imshow(grad(:,col1:min(col1+nbcol,nc)),[])
end