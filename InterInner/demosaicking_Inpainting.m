function X = demosaicking_Inpainting(img)

green = GetBayerForTV(img);
Xc = green;

load('mask_g');
% load('mask_g_imax');
% load('mask_g_t');
% load('mask_g_test');

B = Xc;
B(B<0) = 0; B(B>255) = 255; % Pixels in range 0,...,255.

mask=mask_g;
% mask=mask_g_imax;
% mask = mask_g_t;
% mask = mask_g_test;

mask(mask>0) = 1;                         % Set mask values to 1.
B(mask>0) = 255;                       % Set corrupted pixels to white.

[X,info] = TVinpaint(B,mask,0);
G= Xc + imfilter(Xc, [0 1 0; 1 0 1; 0 1 0]/4);

X = round(0.38*G+0.62*X);
X(X>255)=255;
