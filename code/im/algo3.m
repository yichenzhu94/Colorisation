function [u, E] = algo3(c,sigma,tau,lambda, epsilon)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
   [sx, sy, sd, sc] = size(c);
   
   % Init
   W = ones(sx,sy,sc) ./ sc;
   u = sum(c,4)./ sc;
   ubar = u;
   [pUx, pUy] = grad(u(:,:,1));
   [pVx, pVy] = grad(u(:,:,2));
   cpt = 0;
   
   verbose = true;
   E = [];
   
   % Main loop
   while verbose
       
       % update p
      [gdUx, gdUy] = grad(ubar(:,:,1));
      [gdVx, gdVy] = grad(ubar(:,:,2));
       pUx = pUx + sigma*gdUx;
       pUy = pUy + sigma*gdUy;
       pVx = pVx + sigma*gdVx;
       pVy = pVy + sigma*gdVy;
       
       % projection
       nrm = pUx.*pUx + pUy.*pUy + pVx.*pVx + pVy.*pVy ;
       id = nrm > 1; 
       nrm = sqrt(nrm(id));
       pUx(id) = pUx(id) ./ nrm;
       pUy(id) = pUy(id) ./ nrm;
       pVx(id) = pVx(id) ./ nrm;
       pVy(id) = pVy(id) ./ nrm;
       
       % compute divergence
       dU = div(pUx,pUy);
       dV = div(pVx,pVy);
       
       % devide u into U,V channels
       uU = u(:,:,1);
       uV = u(:,:,2);
       
       % compute the closest candidate for ubar
       diff_old = 10000000;
       for i = 1:sc
           c_vecu = c(:,:,1,i);
           c_vecv = c(:,:,2,i);
           ubar_u = ubar(:,:,1);
           ubar_v = ubar(:,:,2);
           ubar_vec = [ubar_u(:) ubar_v(:)];
           c_vec = [c_vecu(:) c_vecv(:)];
           diff = norm((c_vec - ubar_vec),2);
           if diff < diff_old
               C_idx = i;
           end
           diff_old = diff;
       end
       CU = c(:,:,1,C_idx);
       CV = c(:,:,2,C_idx);
       % update u
       uU = (uU + tau*(dU + lambda*CU))/(1+tau*lambda);
       uV = (uV + tau*(dV + lambda*CV))/(1+tau*lambda);
       
%        uU(uU < mU ) = mU;
%        uU(uU > MU ) = MU;
%        uV(uV < mV ) = mV;
%        uV(uV > MV ) = MV;
       
       % Stop criterier
       uUold = u(:,:,1);
       uVold = u(:,:,2);
       nrm = (uUold -uU) .* (uUold - uU) + (uVold - uV) .* (uVold-uV) ;
       nrm = sum(sum(nrm));
       if nrm < epsilon 
           verbose = false ;
       end
       display(nrm);
       
       unew(:,:,1)  = uU;
       unew(:,:,2)  = uV;
       % Update ubar
       ubar = 2*unew - u;
       u = unew;
       
       % Energy
       [DUx,DUy] = grad(uU);
       [DVx,DVy] = grad(uV);
       tv = sum(sum(sqrt(DUx.^2 + DUy.^2 + DVx.^2 + DVy.^2)));
       E = [E,tv];
       cpt = cpt + 1;
       fprintf("iteration %d : E = %.10g\n",cpt,E(end));
   end

end

