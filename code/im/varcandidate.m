function V = varcandidate(S,Sm, sigma, coef)
    %computes the varpatch value
    %for each candidate pixel in the subsampled source image
    %S = source image
    %Sm = subsampled sourc image
    %Sigma = size of the window of the varpatch
    %coef = step fo the subsampling
    
    [nx, ny] = size(Sm);
    [nr, nc] = size(S);
    V = zeros(nx,ny);
    %V = zeros(sigma+1,sigma+1,nx, ny);
    
    for rm = 1:nx
        r = (rm-1)*coef+1;
        if ( (r > round(sigma/2)) & (r < nr-round(sigma/2) ))
            for cm =1:ny
                c = (cm-1)*coef+1;
                if ( (c > round(sigma/2)) & (c < nc-round(sigma/2)) )
                    V(rm,cm) = varpatch(r,c,sigma, S);
                end                
            end

        end
    end
    
    
end

