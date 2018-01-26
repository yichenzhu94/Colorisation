function [Dx Dy] = grad(u)
%taken from__ Inverse problem with 
%quadratic data-fidelity and total variation regularity__Remi Arbegel

[ny, nx] = size(u);
 Dx = u(:,[2:nx,nx]) - u;
 Dy = u([2:ny,ny],:) - u;
end
