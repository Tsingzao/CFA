function green = GetBayerForTV(img)

img = double(img);

load('mask_g');
% load('mask_g_imax');
% load('mask_g_t');
% load('mask_g_test');

green = img(:,:,2);

green = green .* (1 - mask_g);
% green=green.*(1-mask_g_imax);
% green = green .* (1 - mask_g_t);
% green=green.*(1-mask_g_test);