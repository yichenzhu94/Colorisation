function [qx, qy] = varmatch(T, Sc, coef, x, y, sigma)
% T is the image from wigh we take the pixel value
%(x,y) are the coordinates of the point int T
%Sc are the candidates variance
%coef of subsampling ( so we have the mapping between Sc coordinates and
%the initial coordinates)
%sigma is the size of the patch

[nx, ny] = size(Sc);
%[nr, nc] = size(S);
%We first estimate the target variance
vobj = varpatch(x,y,sigma, T);
% min = 10*vobj;
% qx=0; qy = 0;
Sc = abs(Sc-vobj);
[qx, qy] = find(Sc == min(Sc(:)));
l = length(qx);
idx = randsample(1:l,1);
qx = coef*(qx(idx,1)-1) + 1;
qy = coef*(qy(idx,1)-1) + 1;

% for each candidate in Sm we compare the variance(vtmp) to the the target
% variance and we retain the coordinate ( in S)  of the point that has the
% closest variance to the target value.

% for rm = 1:nx
%     r = (rm-1)*20+1;
%     if ( (r > round(sigma/2)) & (r < nr-round(sigma/2) ))
%         for cm =1:ny
%             c = (cm-1)*20+1;
%             if ( (c > round(sigma/2)) & (c < nc-round(sigma/2)) )
%                 vtmp = varpatch(r,c,sigma, S);
%                 if ( abs(vobj-vtmp) < min)
%                     min = abs(vobj-vtmp);
%                    qx =r; qy = c; 
%                 end                
%             end
% 
%         end
%     end
%     
%     
% end

end