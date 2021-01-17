function [Infl1,Infl2,Gradi1,Gradi2]=inflexImg (bScan_p)
% INFLEXIMG find the weights of the inflexion points in the image BSCAN_P
% BSCAN_P : is one Bscan
%   INFL1 : weight of the inflexion points where the gradient is positive(weight = 0 ==> no inflexion
%   point)
%   INFL2 : weight of the inflexion points where the gradient is
%   negative
%   GRAD1 : positive gradient
%   GRAD2 : negative gradient
    filteredBscan = single(mat2gray(bScan_p));%normalisation of the image (might be useless)

    tg = 2;
    mask_grad = (heaviside([-tg:+tg])-0.5)';%mask of the gradient
   
    th_grad =0;
   
    Gradi = imfilter(filteredBscan,mask_grad,'symmetric');%calculus of the gradient

    grad2 = 0.25*imfilter(filteredBscan,[1;-2;1],'symmetric');%calculus of the laplacian
    lap2 = [zeros(1,size(grad2,2),'single');grad2(1:end-1,:)];
    mult = grad2.*lap2;
    signC = mult<0; %find where the laplacian sign change (=definition of the inflexion point)
    
    Infl1 =  (Gradi > th_grad) & signC; %position of the positive inflexion points
    Infl2 =  (Gradi < th_grad) & signC; %position of the negative inflexion points

    Infl1 = Infl1.*Gradi; %we add the weight
    Infl2 = abs(Infl2.*Gradi);%we add the weight
    Gradi1 = Gradi; Gradi1(Gradi1<0)=0;%dark to bright
    Gradi2 = -Gradi; Gradi2(Gradi2<0)=0;%bright to dark
Infl1(end-2:end,:) = 0;Infl2(end-2:end,:) = 0;%delete the first and last lines where there are unwanted border effects
Gradi1(end-2:end,:) = 0;Gradi2(end-2:end,:) = 0;

%figure; plot(Gradi(:,251))
% figure;subplot(3,1,1),imshow(filteredBscan,[]), subplot(3,1,2),imshow(Gradi1,[]),axis on,grid on, subplot(3,1,3),imshow(Gradi2,[]),axis on,grid on;
% figure; subplot(2,1,1),imshow(Infl1,[]),axis on,grid on, subplot(2,1,2),imshow(Infl2,[]),axis on,grid on;

end

