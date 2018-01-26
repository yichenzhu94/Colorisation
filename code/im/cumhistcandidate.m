function Hist = cumhistcandidate(S,Sm, sigma, coef)
    %computes the cumhistpatch value
    %for each candidate pixel in the subsampled source image
    %S = source image
    %Sm = subsampled sourc image
    %Sigma = size of the window of the varpatch
    %coef = step fo the subsampling
    [nx, ny] = size(Sm);
    [nr, nc] = size(S);
    Hist = zeros(255,nx, ny);
    
    for rm = 1:nx
        r = (rm-1)*coef+1;
        if ( (r > round(sigma/2)) && (r < nr-round(sigma/2) ))
            for cm =1:ny
                c = (cm-1)*coef+1;
                if ( (c > round(sigma/2)) && (c < nc-round(sigma/2)) )
                    
                 
                    Hist(:,rm,cm) = cumhistpatch(r,c,sigma, S);
                end
            end

        end
    end
    
end