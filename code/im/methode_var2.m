%           Primal-dual algorithm for coloriztion with TV regularization
%                       
%  
%   c = candidates color : dim = (height,width,2, numcandidate)
%   lambda = regularity parameter (proximity of u to a candidate)
%   alpha = regularity term ( sparsity of coef w)
%   epsilon = stopping criterion
%   tau = learning steps for z, x and w
%   u = restored image
%   E = energy (E(i) = energy computed at iteration i)
%  

function [u, E] = methode_var2(c,lambda,sigma,rho,tau, Y, epsilon)


%     retrieve dimensions of domains omega & Omega (recall
%     that operator A is from R^Omega to R^omega).
  
   s = size(c);
   sx = s(1);
   sy = s(2);
   sd = s(3);
   sc = s(4);
   mU = -0.436 * 255;
   mV = -0.615 * 255;
   MU = 0.436 * 255;
   MV = 0.615 * 255;

   % initialize primal and dual variables
   W = ones(sx,sy,sc) ./ sc;
   u = sum(c,4)./ sc;
   ubar = u;
   z = zeros(sx,sy,4);
   zUx = z(:,:,1) ;
   zUy = z(:,:,2);
   zVx = z(:,:,3) ;
   zVy = z(:,:,4);
   ctr = 0;
   
   iter = true ;
   E = [];
   [DYx,DYy] = grad(Y);
   
   %% **************************** main loop ****************************//
   while iter
     % compute ((DUx,DUy),(DVx,DVy) ) = (grad(uU), grad(uV)
     uU = ubar(:,:,1);
     uV = ubar(:,:,2);
%      uU = u(:,:,1);
%      uV = u(:,:,2);
     [DUx,DUy] = grad(uU);
     [DVx,DVy] = grad(uV);
     % update z = (zUx,zUy, zVx, zVy) 
     zUx  = zUx + sigma*DUx;
     zUy  = zUy + sigma*DUy;
     zVx  = zVx + sigma*DVx;
     zVy  = zVy + sigma*DVy;
     nrm = zUx.*zUx + zUy.*zUy + zVx.*zVx + zVy.*zVy + DYx .* DYx + DYy .* DYy ;
     id = nrm > 1; nrm = sqrt(nrm(id));
     zUx(id) = zUx(id) ./ nrm;
     zUy(id) = zUy(id) ./ nrm;
     zVx(id) = zVx(id) ./ nrm;
     zVy(id) = zVy(id) ./ nrm;
     
    
     
     % update W
     nrm = zeros(sx,sy,sc);
     for i = 1:sc
         tmp = (ubar-c(:,:,:,i)) .* (ubar-c(:,:,:,i));
         tmp = sum(tmp,3);
         
         nrm(:,:,i) = tmp;
     end
     nrm = rho*lambda*nrm  ;
     W = (W - nrm);
     for i = 1:sx
        for j = 1:sy
           W(i,j,:) = projsplx( W(i,j,:) ) ;
        end
     end
     %projection into the canonical basis 
     for i = 1:sx
        for j = 1:sy
            tmp = W(i,j,:);
            [~, idx] = sort( tmp(:) ) ;
            tmp = zeros(1,1,sc);
            tmp(1,1,idx(end)) = 1;
            W(i,j,:) = tmp;
            
            
        end
     end
     
      % compute d = div(p) 
     dU = div(zUx,zUy);
     dV = div(zVx,zVy);
     
     % update u 
     uU = u(:,:,1);
     uV = u(:,:,2);
%      [DUx,DUy] = grad(uU);
%      [DVx,DVy] = grad(uV);
     ctmp = permute(c,[1,2,4,3]);
     cU = ctmp(:,:,:,1);
     cV = ctmp(:,:,:,2);
     tmpU = cU .* W;
     tmpU = sum(tmpU,3);
     tmpV = cV .* W;
     tmpV = sum(tmpV,3);
     uU =  (uU + tau *( dU + lambda * tmpU ))/(1+tau*lambda);
     uV =  (uV + tau *(dV + lambda * tmpV ))/(1+tau*lambda);
     
     %we make sure that the values of U ans V are authorized
     uU(uU < mU ) = mU;
     uU(uU > MU ) = MU;
     uV(uV < mV ) = mV;
     uV(uV > MV ) = MV;
     
     % tcheck if it is necessary to have another iteration
     uUold = u(:,:,1);
     uVold = u(:,:,2);
     norm = (uUold -uU ) .* (uUold - uU) + (uVold - uV) .* (uVold-uV) ;
     norm = sum(sum(norm));
     if norm < epsilon 
         iter = false ;
     end
     display(norm);
     
     %update ubar
     ubarU = (2*uU - uUold );
     ubarV = (2*uV - uVold );
     
    
     u(:,:,1)  = uU;
     u(:,:,2)  = uV;
     ubar(:,:,1)  = ubarU;
     ubar(:,:,2)  = ubarV;
     % compute energy of u
     [DUx,DUy] = grad(uU);
     [DVx,DVy] = grad(uV);
     tv = sum(sum(sqrt(DUx.^2 + DUy.^2 + DVx.^2 + DVy.^2 + DYx .^2 + DYy.^2 )));
%      sparse = (-1*W + 1) .* W;
     %sparse = 0.5 * alpha * sum(sum(sum(sparse,3)));
     
     nrm = zeros(sx,sy,sc);
     for i = 1:sc
         tmp = (u-c(:,:,:,i)) .* (u-c(:,:,:,i));
         tmp = sum(tmp,3);
         nrm(:,:,i) = tmp;
     end
     nrm = nrm .* W ;
     dataterm = 0.5 * lambda * sum(sum(sum(nrm,3)));
     
     E = [E, dataterm + tv];
     ctr = ctr + 1;
     
     fprintf("iteration %d : E = %.10g\n",ctr,E(end));
   end
   %% ************************* end of main loop *************************//
 end