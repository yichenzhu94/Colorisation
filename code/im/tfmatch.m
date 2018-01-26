function [qx, qy] = tfmatch(T, Sc, coef, x, y, sigma)
[~, ~, nx, ny] = size(Sc);
%[nr, nc] = size(S);
%We first estimate the target variance
tfobj = tfpatch(x,y,sigma, T);
% min = 10*vobj;
% qx=0; qy = 0;
tmp = zeros(nx,ny);
% display(size(tfobj));
% display(size(Sc));
Sc = sum(sum(abs(Sc-tfobj)));
tmp(:,:)=Sc(1,1,:,:);
[qx, qy] = find(tmp == min(tmp(:)));
qx = coef*(qx(1,1)-1) + 1;
qy = coef*(qy(1,1)-1) + 1;

end