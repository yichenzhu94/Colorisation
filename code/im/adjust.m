function Tm = adjust(S,T)
%apply affine transformation to the target image so it has same variance
%and mean as the source function
%both images are in grey scale
ms = mean(S(:));
ss = std(S(:));
mt = mean(T(:));
st = std(T(:));
Tm = (ss/st)*T + ms - (ss/st)*mt;
end