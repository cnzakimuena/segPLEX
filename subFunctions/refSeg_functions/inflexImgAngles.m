function [infl_out1,infl_out2,grad_out1,grad_out2] = inflexImgAngles(volum,angles)
%%INFLEXIMGANGLES calculates the positions and the 'weights' of the inflexion points in an image.
%Here, the inflexion points are the points where the sign of the second derivative
%(according to a given angle) changes.
%What we call 'weight of inflexion points' is the value of the gradient(first derivative) at
%the inflexion point. 
%For each image, we calculate the gradient according to an angle. angle = 0 is
%the gradient on each A-scan. angle = 90 is the gradient on the line
%If the variable ANGLES contains several angles, the weights of the
%inflexion points are the mean of the weight for each angle.
%There is two types of inflexion points : those whose gradient is positive
%(INFL_OUT1) and those whose gradient is negative (INFL_OUT2)
%
%   VOLUM : 3D OCT volume
%   ANGLES : value of the angles on which we calculate the gradient
%   INFL_OUT1 : weight of the inflexion points where the gradient is positive(weight = 0 ==> no inflexion
%   point)
%   INFL_OUT1 : weight of the inflexion points where the gradient is
%   negative
%   GRAD_OUT1 : positive gradient
%   GRAD_OUT2 : negative gradient

[nbL,nbC,nbZ] = size(volum);
nbPadv = 30;
nbPadh = round(nbC/10);
volum2= padarray(volum,[nbPadv, nbPadh,0],'symmetric');
nbAngles  = length(angles);
[nbL,nbC,nbZ] = size(volum2);
infl_out2_1 = zeros(nbL,nbC,nbZ,'single');
infl_out2_2 = zeros(nbL,nbC,nbZ,'single');
grad_out2_1 = zeros(nbL,nbC,nbZ,'single');
grad_out2_2 = zeros(nbL,nbC,nbZ,'single');
%infl_out2_2 = zeros(nbL,nbC,nbZ);
%figure; imshow(img,[])

for ind_bScan = 1:nbZ %for each bscan
    cpt = 1;
    img = volum2(:,:,ind_bScan);
    
    infl_out1 = zeros(size(img,1),size(img,2),2*nbAngles-1,'single');
    infl_out2 = zeros(size(img,1),size(img,2),2*nbAngles-1,'single');
    grad_out1 = zeros(size(img,1),size(img,2),2*nbAngles-1,'single');
    grad_out2 = zeros(size(img,1),size(img,2),2*nbAngles-1,'single');
    
    for i = 1:nbAngles %for each angles, find inflexion points with there weights
        
        rotation  = imrotate(img,angles(i),'bilinear');%rotation of the image according to the angle(if angle =0 : useless)
        
        rotation_ = imrotate(img,-angles(i),'bilinear');
    
        [Infl1,Infl2,Gradi1,Gradi2]      = inflexImg (rotation); 
        
        [Infl_1,Infl_2,Gradi_1,Gradi_2]     = inflexImg (rotation_);
        %figure; subplot(2,1,1),imshow(logical(Infl),[]),subplot(2,1,2),imshow(Infl,[])
        %figure; subplot(2,1,1),imshow(logical( Infl_ ),[]),subplot(2,1,2),imshow( Infl_ ,[])
        
        %rotation of the image according in the opposite angle (if the
        %angle =0) useless
        resultats  = antirotation(cat(3,Infl1,Infl2,Gradi1,Gradi2) ,angles(i) ,[size(img)],1);
       
        resultats_  = antirotation(cat(3,Infl_1,Infl_2,Gradi_1,Gradi_2) ,-angles(i) ,[size(img)],1);
     
        %figure; subplot(3,1,1),imshow(img,[]),subplot(3,1,2),imshow(Infl,[]),subplot(3,1,3),imshow(Infl_,[])
        %figure; imshow(rotation_one)
        %figure,subplot(4,1,1),imshow(resultats1,[]),subplot(4,1,2),imshow(resultats2(2:end,:),[]),subplot(4,1,3),imshow(resultats_1(2:end,:),[]),subplot(4,1,4),imshow(resultats_2(2:end,:),[])
        
        %the rotation might have modified the size of the image. A resizing
        %is needed
        resultats  = imresize(resultats,[size(img,1),size(img,2)],'nearest');
        resultats_  = imresize(resultats_,[size(img,1),size(img,2)],'nearest');

        if mod(angles(i),360) ==0
            infl_out1(:,:,cpt) = resultats(:,:,1);
            infl_out2(:,:,cpt) = resultats(:,:,2);
            grad_out1(:,:,cpt) = resultats(:,:,3);
            grad_out2(:,:,cpt) = resultats(:,:,4);
            cpt = cpt+1;
        else
            infl_out1(:,:,cpt) = resultats(:,:,1);
            infl_out1(:,:,cpt+1) = resultats_(:,:,1);
            grad_out1(:,:,cpt) = resultats(:,:,3);
            grad_out1(:,:,cpt) = resultats_(:,:,3);
            
            infl_out2(:,:,cpt) = resultats(:,:,2);
            infl_out2(:,:,cpt+1) = resultats_(:,:,2);
            grad_out2(:,:,cpt) = resultats(:,:,4);
            grad_out2(:,:,cpt) = resultats_(:,:,4);
            cpt = cpt+2;
        end
      
        
    end
    
    %calculus of the mean of the weights for each angles
    infl_out1 = mean(infl_out1,3);
    infl_out1(1:2,:) = 0;
    infl_out2 = mean(infl_out2,3);
    infl_out2(1:2,:) = 0;
    
    grad_out1 = mean(grad_out1,3);
    grad_out1(1:2,:) = 0;
    grad_out2 = mean(grad_out2,3);
    grad_out2(1:2,:) = 0;
    
    
    infl_out2_1(:,:,ind_bScan) = infl_out1;
    infl_out2_2(:,:,ind_bScan) = infl_out2;
    
    grad_out2_1(:,:,ind_bScan) = grad_out1;
    grad_out2_2(:,:,ind_bScan) = grad_out2;
    %figure; subplot(2,1,1),imshow(img_out1,[]),subplot(2,1,2),imshow(img_out2,[])
end
infl_out1=infl_out2_1(nbPadv+1:end-nbPadv,nbPadh+1:end-nbPadh,:);
infl_out2=infl_out2_2(nbPadv+1:end-nbPadv,nbPadh+1:end-nbPadh,:);

grad_out1=grad_out2_1(nbPadv+1:end-nbPadv,nbPadh+1:end-nbPadh,:);
grad_out2=grad_out2_2(nbPadv+1:end-nbPadv,nbPadh+1:end-nbPadh,:);
%figure; imshow3D(grad_out2,[])
grad_out1(end-1:end,:,:) = 0;
grad_out2(end-1:end,:,:) = 0;

end

