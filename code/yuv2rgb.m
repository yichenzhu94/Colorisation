function rgb = yuv2rgbb(YUV)
%takes an image yuv and returns an image in the rgb colorspace

% invD is the matrix correspond to the linear transformation yuv -> rgb
invD = [1 0 1.13983;
        1 -0.39465 -0.58060;
        1 2.03211 0];

% apply the transformation invD on each channel
y = YUV(:,:,1);
u = YUV(:,:,2);
v = YUV(:,:,3);

[nr nc ~] = size(YUV);
YUV = [y(:)';u(:)';v(:)'];
RGB = invD*YUV;


% we express the image in the initial dimention
r = RGB(1,:);
g = RGB(2,:);
b = RGB(3,:);
r = reshape(r', [nr nc]);
g = reshape(g', [nr nc]);
b = reshape(b', [nr nc]);
rgb = zeros(nr, nc, 3);
rgb(:,:,1) = r; rgb(:,:,2) = g ; rgb(:,:,3) = b;
end