function GetImgMask(img)

tic;
fprintf( sprintf('\nRunning demosaick CombinedAlgorithm ... \n'));
rgb = double(img);
pattern = 'grbg';

%% mosaic and mask
[mosaic mask] = mosaic_bayer(rgb, pattern);
%%
A = extendImage(rgb,10);

[M,N,ch]=size(A);
h = zeros(M,N,3);
g = zeros(M,N);
if ch==3
    h(1:2:M,1:2:N,2)=A(1:2:M,1:2:N,2);
    h(2:2:M,2:2:N,2)=A(2:2:M,2:2:N,2);
    h(1:2:M,2:2:N,1)=A(1:2:M,2:2:N,1);
    h(2:2:M,1:2:N,3)=A(2:2:M,1:2:N,3);
    g(1:2:M,1:2:N)=A(1:2:M,1:2:N,2);
    g(2:2:M,2:2:N)=A(2:2:M,2:2:N,2);
    g(1:2:M,2:2:N)=A(1:2:M,2:2:N,1);
    g(2:2:M,1:2:N)=A(2:2:M,1:2:N,3);
else
    g = A;
    h(1:2:M,1:2:N,2)=A(1:2:M,1:2:N);
    h(2:2:M,2:2:N,2)=A(2:2:M,2:2:N);
    h(1:2:M,2:2:N,1)=A(1:2:M,2:2:N);
    h(2:2:M,1:2:N,3)=A(2:2:M,1:2:N);
end

f=[-1/4 1/2 1/2 1/2 -1/4];
Ah=conv2(g,f);
Ah=Ah(:,3:2+N);
Av=conv2(g,f');
Av=Av(3:2+M,:);

dh=zeros(M,N);
dv=dh;
for i=1:2:M
   dh(i,1:2:N)=g(i,1:2:N)-Ah(i,1:2:N);
   dh(i,2:2:N)=Ah(i,2:2:N)-g(i,2:2:N);
   dv(i,1:2:N)=g(i,1:2:N)-Av(i,1:2:N);
   dv(i,2:2:N)=Av(i,2:2:N)-g(i,2:2:N);
end
for i=2:2:M
   dh(i,2:2:N)=g(i,2:2:N)-Ah(i,2:2:N);
   dh(i,1:2:N)=Ah(i,1:2:N)-g(i,1:2:N);
   dv(i,2:2:N)=g(i,2:2:N)-Av(i,2:2:N);
   dv(i,1:2:N)=Av(i,1:2:N)-g(i,1:2:N);
end

vmap = zeros(M,N);
hmap = zeros(M,N);

for i=4:1:M-4
    for j=4:1:N-4
        vmap(i,j) = abs(dv(i-1,j)-dv(i+1,j))+abs(dv(i-1,j-1)-dv(i+1,j-1))+abs(dv(i-1,j+1)-dv(i+1,j+1));
        hmap(i,j) = abs(dh(i,j-1)-dh(i,j+1))+abs(dh(i-1,j-1)-dh(i-1,j+1))+abs(dh(i+1,j-1)-dh(i+1,j+1));
    end
end

vmapu = zeros(M,N);
hmapl = zeros(M,N);
vmapd = zeros(M,N);
hmapr = zeros(M,N);

vmapv = zeros(M,N);
hmaph = zeros(M,N);

r55 = fspecial('gaussian',5,2);  
for i=5:1:M-5
    for j=5:1:N-5    
        vmapu(i,j) = 100/(sum(sum(vmap(i-4:i,j-2:j+2).*r55))^2+1);
        vmapd(i,j) = 100/(sum(sum(vmap(i:i+4,j-2:j+2).*r55))^2+1);
        hmapl(i,j) = 100/(sum(sum(hmap(i-2:i+2,j-4:j).*r55))^2+1);
        hmapr(i,j) = 100/(sum(sum(hmap(i-2:i+2,j:j+4).*r55))^2+1);
        
        vmapv(i,j) = 100/(sum(sum(vmap(i-2:i+2,j-2:j+2).*r55))^2+1);
        hmaph(i,j) = 100/(sum(sum(hmap(i-2:i+2,j-2:j+2).*r55))^2+1);
    end
end

x1 = [12 17 21 24 26]/100;
x2 = [26 24 21 17 12]/100;

for i=6:2:M-6
    for j=6:2:N-6
        h(i+1,j,2) = ((sum(dv(i-3:i+1,j).*x1'))*vmapu(i+1,j) + (sum(dv(i+1:i+5,j).*x2'))*vmapd(i+1,j) + (sum(dh(i+1,j-4:j).*x1))*hmapl(i+1,j) + (sum(dh(i+1,j:j+4).*x2))*hmapr(i+1,j))/(vmapu(i+1,j)+vmapd(i+1,j)+hmapl(i+1,j)+hmapr(i+1,j)) + g(i+1,j);
    end
end

for i=6:2:M-6
    for j=6:2:N-6
        h(i,j+1,2) = ((sum(dv(i-4:i,j+1).*x1'))*vmapu(i,j+1) + (sum(dv(i:i+4,j+1).*x2'))*vmapd(i,j+1) + (sum(dh(i,j-3:j+1).*x1))*hmapl(i,j+1) + (sum(dh(i,j+1:j+5).*x2))*hmapr(i,j+1))/(vmapu(i,j+1)+vmapd(i,j+1)+hmapl(i,j+1)+hmapr(i,j+1)) + g(i,j+1);
    end
end
green = h(11:M-10,11:N-10,2);
%%
% Red and Blue demosaicking
%%
% red = red_interpolation(green, mosaic, mask, 5, 5, 0);
for i=4:2:M-4
    for j=4:2:N-4
        h(i,j+1,1) = h(i,j+1,2) - ((h(i-1,j,2)-h(i-1,j,1))+(h(i-1,j+2,2)-h(i-1,j+2,1))+(h(i+1,j,2)-h(i+1,j,1))+(h(i+1,j+2,2)-h(i+1,j+2,1)))*1.25/4 ...
            + ((h(i-3,j,2)-h(i-3,j,1))+(h(i-1,j-2,2)-h(i-1,j-2,1))+(h(i-3,j+2,2)-h(i-3,j+2,1))+(h(i-1,j+4,2)-h(i-1,j+4,1))+(h(i+3,j,2)-h(i+3,j,1))+(h(i+1,j-2,2)-h(i+1,j-2,1))+(h(i+3,j+2,2)-h(i+3,j+2,1))+(h(i+1,j+4,2)-h(i+1,j+4,1)))*0.25/8;
    end
end
for i=6:2:M-6
    for j=6:2:N-6
        h(i,j,1) = h(i,j,2) - (((h(i,j-1,2)-h(i,j-1,1)) + (h(i,j+1,2)-h(i,j+1,1)))*hmaph(i,j) + ((h(i-1,j,2)-h(i-1,j,1)) + (h(i+1,j,2)-h(i+1,j,1)))*vmapv(i,j))/(2*(vmapv(i,j) + hmaph(i,j)));
        h(i+1,j+1,1) = h(i+1,j+1,2) - (((h(i,j+1,2)-h(i,j+1,1)) + (h(i+2,j+1,2)-h(i+2,j+1,1)))*vmapv(i+1,j+1) + ((h(i+1,j,2)-h(i+1,j,1)) + (h(i+1,j+2,2)-h(i+1,j+2,1)))*hmaph(i+1,j+1))/(2*(vmapv(i+1,j+1) + hmaph(i+1,j+1)));
    end
end
red = h(11:M-10,11:N-10,1);
%%
blue  = blue_interpolation(green, mosaic, mask, 5, 5, 0);
%%
% result image
rgb_dem(:,:,1) = red;
rgb_dem(:,:,2) = green;
rgb_dem(:,:,3) = blue;
%%
% calculate PSNR and CPSNR
psnr = impsnr(rgb, rgb_dem, 255, 10);
cpsnr = imcpsnr(rgb, rgb_dem, 255, 10);

% print PSNR and CPSNR
disp('__________________________________')
fprintf( sprintf( 'Red:%f\n',     psnr(1)     ) );
fprintf( sprintf( 'Green:%f\n',   psnr(2)     ) );
fprintf( sprintf( 'Blue:%f\n',    psnr(3)     ) );
fprintf( sprintf( '::::: CPSNR ::::::::%f\n',   cpsnr       ) );

toc;