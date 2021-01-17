function [vol_out,trans_out]=flattenVol(vol,lRPE,ref,err)
%%FLATTENVOL flattens each bscan so that the BM is straight on the ref th
%%line. 

[nbl,nbc,nbz] =  size(vol);
vol_out = single(zeros(size(vol)));

trans_out = zeros(nbz,nbc);
vol_out(:,:,err) = nan;
trans_out(err,:) = nan;
z = 1:nbz;
%setdiff(z,err) returns the data in z that is not in err, 
%with no repetitions
z=setdiff(z,err);
for c = 1:nbc
    for z2=z 
        trans = round(ref-lRPE(z2,c));
        if trans<=0
        vol_out(1:nbl+trans,c,z2) = vol(-trans+1:nbl,c,z2);
        else
        vol_out(trans+1:nbl,c,z2) = vol(1:nbl-trans,c,z2);
        end
        trans_out(z2,c) = trans;
    end
end
end