%% load and rescale images so that they have same dimention
addpath('../../YUV_tools/YUV');
addpath('../../projsplx');
RGB = imread('test.jpg');
%RGB =imresize(RGB, [300 500] );
figure(1);imshow(RGB, []);
display(size(RGB));
RGB2 = imread('test.jpg');
figure(2); imshow(RGB2,[]);
[nr nc ~] = size(RGB);
RGB2 =imresize(RGB2, [nr nc] );
display(size(RGB2));

%% Convert images in YUV
R1 = RGB(:,:,1);
G1 = RGB(:,:,2);
B1 = RGB(:,:,3);
R2 = RGB2(:,:,1);
G2 = RGB2(:,:,2);
B2 = RGB2(:,:,3);
[Y1,U1,V1] = rgb2yuv(R1,G1,B1);
[Y2, U2,V2] = rgb2yuv(R2,G2,B2);
S = double(Y1);
T = double(Y2);

%% adjust the greyscale of the target image so it coresspond 
% to the source image
%Tm = adjust(S, T);
Tm = T;


%% subsampling in the source image
% we dont want to compare every target pixel to every source pixel
% instead, we compare each target to a subsample of the source
%coef is the coefficient of subsampling
coef = 50;
Ss = S(1:coef:end, 1:coef:end);


%% We get the candidate color for each pixel of the taget image
% sigm(0) and sigma(1) = windows for the variance candidates
% sigm(2) sigma(3) and sigma(4) are the windows size for the TF
% sigm(5) sigm(6) and sigm(7) are the windows for the cumulativ histograms

sigm = [6, 8, 6, 8, 10, 6, 8, 10];
[candidate, resultat] = candidat(S, Ss, T,Tm, Y1, U1, V1, sigm, coef);
% candidate = candidate(:,:,:,3:5);
 test1 = resultat(:,:,:,1);
% test2 = resultat(:,:,:,2);
% test3 = resultat(:,:,:,3);
% test4 = resultat(:,:,:,4);
% test5 = resultat(:,:,:,5);
% test6 = resultat(:,:,:,6);
test7 = resultat(:,:,:,7);
% test8 = resultat(:,:,:,8);
% Final = yuv2rgb(test1(:,:,1),test1(:,:,2), test1(:,:,3));
% figure(1);imshow(Final, []);
% Final = yuv2rgb(test2(:,:,1),test2(:,:,2),test2(:,:,3));
% figure(3);imshow(Final, []);
% Final = yuv2rgb(test3(:,:,1),test3(:,:,2),test3(:,:,3));
% figure(3);imshow(Final, []);
% Final = yuv2rgb(test1(:,:,1),test4(:,:,2), test4(:,:,3));
% figure(4);imshow(Final, []);
% Final = yuv2rgb(test2(:,:,1),test5(:,:,2),test5(:,:,3));
% figure(5);imshow(Final, []);
% Final = yuv2rgb(test3(:,:,1),test6(:,:,2),test6(:,:,3));
% figure(6);imshow(Final, []);
Final = yuv2rgb(test1(:,:,1),test7(:,:,2), test7(:,:,3));
figure(3);imshow(Final, []);
% Final = yuv2rgb(test2(:,:,1),test8(:,:,2),test8(:,:,3));
% figure(8);imshow(Final, []);


%% Find the best choice of candidates that minimizes the objectif function
%first method is the classic primal-dual with no coupling
lambda = 1;
alpha = 1;
tauw = 0.45;
tauz = 0.225;
tauu = 0.5;
epsilon = 5;

[u, E, W] = methode_var_origin(candidate,lambda,alpha,tauw,tauz,tauu, epsilon);

u = yuv2rgb(resultat(:,:,1,1),u(:,:,1),u(:,:,2));
figure(6);imshow(u,[]); % display output (restored)
figure(7); plot(1:length(E),E); % plot energy evolution
xlabel("iteration"); ylabel("energy");
title("Energy decrease");
figure(8); hist(W(:));


%second method is a chambolle and Pock algorithm with no coupling and we change
%the projection of the coef w

sigma = 0.005;
tau = 5;
lambda = 0.005;
rho = 1000;
epsilon = 1800;

[u, E, W] = methode_var_nocoupling(candidate,lambda,sigma,rho,tau, epsilon);

u = yuv2rgb(resultat(:,:,1,1),u(:,:,1),u(:,:,2));
figure(9);imshow(u,[]); % display output (restored)
figure(10); plot(1:length(E),E); % plot energy evolution
xlabel("iteration"); ylabel("energy");
title("Energy decrease");
figure(11); hist(W(:));

% 
% %third method : we ignore the coefficient w
% 

[u, E] = algo3(candidate,sigma,tau,lambda, epsilon);

u = yuv2rgb(resultat(:,:,1,1),u(:,:,1),u(:,:,2));
figure(12);imshow(u,[]); % display output (restored)
figure(13); plot(1:length(E),E); % plot energy evolution
xlabel("iteration"); ylabel("energy");
title("Energy decrease");





