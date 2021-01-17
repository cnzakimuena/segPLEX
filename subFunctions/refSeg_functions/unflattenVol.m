
function [volOG, RPEt,RPEb,RVIf] = unflattenVol(vol_f,RPEt_in, RPEb_in, RVI_in,err, trans_out)
%%UNFLATTENVOL unflattens each bscan back to its original appearance

vol_f = single(mat2gray(vol_f));

[nbl,nbc,nbz] =  size(vol_f);
volOG = single(zeros(size(vol_f)));

vol_f(:,:,err) = nan;
trans_out(err,:) = nan;
z = 1:nbz;
z=setdiff(z,err);

% volOG = single(zeros(size(vol_f)));
RPEt = zeros(size(RPEt_in));
RPEb = zeros(size(RPEb_in));
RVIf = zeros(size(RVI_in));

for c = 1:nbc
    for z2=z 
        trans = trans_out(z2,c);
        if trans<=0
            volOG(-trans+1:nbl,c,z2) = vol_f(1:nbl+trans,c,z2);
            RPEt(z2,c) = -trans+RPEt_in(z2,c);
            RPEb(z2,c) = -trans+RPEb_in(z2,c);
            RVIf(z2,c) = -trans+RVI_in(z2,c);
        else    
            volOG(1:nbl-trans,c,z2) = vol_f(trans+1:nbl,c,z2);
            RPEt(z2,c) = trans+RPEt_in(z2,c);
            RPEb(z2,c) = trans+RPEb_in(z2,c);
            RVIf(z2,c) = trans+RVI_in(z2,c);
        end
    end
end
            
% imshow3D(volOG,[])            
% figure;imshow3D(volOG,[],'plot',cat(3,RPEt,RPEb, RVIf),'LineWidth',2) 

end
