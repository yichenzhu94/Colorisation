function [candidate, resultat] = candidat(S, Ss, T, Tm, Y1, U1, V1, sigma, coef)

sigma = sigma(:);
sigmamax = max(sigma);
sigmav = sigma(1:2,1);
sigmaf = sigma(3:5,1);
sigmah = sigma(6:8,1);
[nr nc ~] = size(T);
candidate = zeros(nr-2*sigmamax, nc-2*sigmamax, 2, 8);
resultat = zeros(nr-2*sigmamax, nc-2*sigmamax, 3, 8);
ctr = 0;
for sigm = sigmav'
    ctr = ctr+1;
    Tcropped = T(sigm+1:end-sigm, sigm+1:end-sigm);
    [ctrx, ctry] = size(Tcropped);
    target1 = zeros(ctrx, ctry, 2);
    %target2 = zeros(ctrx, ctry, 2);

    Sc = varcandidate(S,Ss, sigm, coef);
    %Tf = tfcandidate(S,Ss, sigm, coef);
    for i =1:ctrx
        for j = 1:ctry
            %compute the coordinates of the source pixel which is the closest 
            %(in variance or TF) to the target pixel
            x = i+sigm;
            y = j+sigm;
            [q1x, q1y] = varmatch(Tm, Sc, coef, x, y, sigm);
            %[q2x, q2y] = tfmatch(T, Tf, coef, x, y, sigm);
            target1(i,j,1)=q1x; target1(i,j,2)=q1y;
            %target2(i,j,1)=q2x; target2(i,j,2)=q2y;
        end
        %fprintf("row done\n");

    end
     gap = sigmamax - sigm;
%      [tmpx tmpy ~] = size(target1); 
%      target = target1(1+gap: tmpx-gap, 1+gap : tmpy-gap, :);

     
   
    TU1 = zeros(ctrx, ctry);
    TV1 = zeros(ctrx, ctry);

    for i =1:ctrx
        for j = 1:ctry
            x = target1(i,j,1);
            y = target1(i,j,2);

            TU1(i,j)= U1(x,y); TV1(i,j)= V1(x,y);
        end
        
    end
     res = zeros(nr-2*sigmamax, nc-2*sigmamax, 3) ;
     res(:,:,1) = Tcropped(1+gap: ctrx-gap, 1+gap : ctry-gap);
     res(:,:,2) = TU1(1+gap: ctrx-gap, 1+gap : ctry-gap);
     res(:,:,3) = TV1(1+gap: ctrx-gap, 1+gap : ctry-gap);
     size(res)
     resultat(:,:,:,ctr) = res;
     candidate(:,:,:,ctr) = res(:,:,2:3);
     size(resultat)


    
end

for sigm = sigmaf'
    ctr = ctr+1;
    Tcropped = T(sigm+1:end-sigm, sigm+1:end-sigm);
    [ctrx, ctry] = size(Tcropped);
    target1 = zeros(ctrx, ctry, 2);
    %target2 = zeros(ctrx, ctry, 2);

    %Sc = varcandidate(S,Ss, sigm, coef);
    Tf = tfcandidate(S,Ss, sigm, coef);
    for i =1:ctrx
        for j = 1:ctry
            %compute the coordinates of the source pixel which is the closest 
            %(in variance or TF) to the target pixel
            x = i+sigm;
            y = j+sigm;
            [q1x, q1y] = tfmatch(T, Tf, coef, x, y, sigm);
            %[q2x, q2y] = tfmatch(T, Tf, coef, x, y, sigm);
            target1(i,j,1)=q1x; target1(i,j,2)=q1y;
            %target2(i,j,1)=q2x; target2(i,j,2)=q2y;
        end
        %fprintf("row done\n");

    end
     gap = sigmamax - sigm
%      [tmpx tmpy ~] = size(target1); 
%      target = target1(1+gap: tmpx-gap, 1+gap : tmpy-gap, :);
% 
%      candidate(:,:,:,ctr) = target;
   
    TU1 = zeros(ctrx, ctry);
    TV1 = zeros(ctrx, ctry);

    for i =1:ctrx
        for j = 1:ctry
            x = target1(i,j,1);
            y = target1(i,j,2);

            TU1(i,j)= U1(x,y); TV1(i,j)= V1(x,y);
        end
        
    end
     res = zeros(nr-2*sigmamax, nc-2*sigmamax, 3) ;
     res(:,:,1) = Tcropped(1+gap: ctrx-gap, 1+gap : ctry-gap);
     res(:,:,2) = TU1(1+gap: ctrx-gap, 1+gap : ctry-gap);
     res(:,:,3) = TV1(1+gap: ctrx-gap, 1+gap : ctry-gap);
     size(res)
     resultat(:,:,:,ctr) = res;
     candidate(:,:,:,ctr) = res(:,:,2:3);
     size(resultat)


    
end

for sigm = sigmah'
    ctr = ctr+1;
    Tcropped = T(sigm+1:end-sigm, sigm+1:end-sigm);
    [ctrx, ctry] = size(Tcropped);
    target1 = zeros(ctrx, ctry, 2);
    %target2 = zeros(ctrx, ctry, 2);

    %Sc = varcandidate(S,Ss, sigm, coef);
    Th = cumhistcandidate(S,Ss, sigm, coef);
    for i =1:ctrx
        for j = 1:ctry
            %compute the coordinates of the source pixel which is the closest 
            %(in variance or TF) to the target pixel
            x = i+sigm;
            y = j+sigm;
            [q1x, q1y] = cumhistmatch(Tm, Th, coef, x, y, sigm);
            %[q2x, q2y] = tfmatch(T, Tf, coef, x, y, sigm);
            target1(i,j,1)=q1x; target1(i,j,2)=q1y;
            %target2(i,j,1)=q2x; target2(i,j,2)=q2y;
        end
        %fprintf("row done\n");

    end
     gap = sigmamax - sigm;
%      [tmpx tmpy ~] = size(target1); 
%      target = target1(1+gap: tmpx-gap, 1+gap : tmpy-gap, :);

     
   
    TU1 = zeros(ctrx, ctry);
    TV1 = zeros(ctrx, ctry);

    for i =1:ctrx
        for j = 1:ctry
            x = target1(i,j,1);
            y = target1(i,j,2);

            TU1(i,j)= U1(x,y); TV1(i,j)= V1(x,y);
        end
        
    end
     res = zeros(nr-2*sigmamax, nc-2*sigmamax, 3) ;
     res(:,:,1) = Tcropped(1+gap: ctrx-gap, 1+gap : ctry-gap);
     res(:,:,2) = TU1(1+gap: ctrx-gap, 1+gap : ctry-gap);
     res(:,:,3) = TV1(1+gap: ctrx-gap, 1+gap : ctry-gap);
     size(res)
     resultat(:,:,:,ctr) = res;
     size(resultat)
     candidate(:,:,:,ctr) = res(:,:,2:3);


    
end



end