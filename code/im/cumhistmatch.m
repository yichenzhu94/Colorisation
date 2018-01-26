function [qx, qy] = cumhistmatch(T, Sc, coef, x, y, sigma)
[~, nx, ny] = size(Sc);
%[nr, nc] = size(S);
%We first estimate the target variance
chobj = cumhistpatch(x,y,sigma,T);
% min = 10*vobj;
% qx=0; qy = 0;
tmp = zeros(nx,ny);
% display(size(tfobj));
% display(size(Sc));
Sc = sum(abs(Sc-chobj));
tmp(:,:)=Sc(1,:,:);
[qx, qy] = find(tmp == min(tmp(:)));
l = length(qx);
idx = randsample(1:l,1);
qx = coef*(qx(idx,1)-1) + 1;
qy = coef*(qy(idx,1)-1) + 1;

end