function yuv = rgb2yuvb(RGB)
%takes an image rgb and returns an image in the yuv colorspace

% D is the matrix correspond to the linear transformation rgb -> yuv
D = [0.299 0.587 0.114;
    -0.14713 -0.28886 0.436;
    0.615 -0.51498 -0.10001];

% apply the transformation D on each channel
rgb1 = RGB(:,:,1);
rgb2 = RGB(:,:,2);
rgb3 = RGB(:,:,3);
[nr nc ~] = size(RGB);
RGB = [rgb1(:)';rgb2(:)';rgb3(:)'];
YUV = D*double(RGB);

% we express the image in the initial dimention
y = YUV(1,:);
u = YUV(2,:);
v = YUV(3,:);
y = reshape(y', [nr nc]);
u = reshape(u', [nr nc]);
v = reshape(v', [nr nc]);
yuv = zeros(nr, nc, 3);
yuv(:,:,1) = y; yuv(:,:,2) = u ; yuv(:,:,3) = v;
end