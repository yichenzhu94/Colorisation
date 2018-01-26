function chp = cumhistpatch(x,y,sigma, im)
%compute the TF of the patch centered in (x,y) of size sigma in image
%im. im is in grayscale
patch = im( x-(sigma/2): x+(sigma/2), y-(sigma/2): y+(sigma/2));
nbin = (sigma+1)^2;
X = 1:255;
%X = X(2:end);
chp = cumsum(hist(patch(:),X))./nbin;
chp = chp(:);
end
