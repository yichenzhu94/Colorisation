function v = varpatch(x,y,sigma, im)
%compute the variance of the patch centered in (x,y) of size sigma in image
%im. im is in grayscale
patch = im( x-(sigma/2): x+(sigma/2), y-(sigma/2): y+(sigma/2));
v = var(patch(:));
end