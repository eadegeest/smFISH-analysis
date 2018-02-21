function[imnew,imsub] = fftfilter_nowrap_fast_pad(imraw,epsin)
% function[imnew,imsub] = fftfilter_nowrap_b(imraw,epsin)
% 
% INPUT parameters
% imraw : raw image
% epsin : parameter epsilon, all structures on length scales > epsin will
% be removed (= spatial frequencies below 1/epsin)
%
% OUTPUT parameters
%
% imnew : the background image (= only structures > epsin)
% imsub : raw image - background image


imraw = double(imraw);
Nraw = size(imraw);

Npad = 50;
tempim = zeros(size(imraw,1)+2*Npad,size(imraw,2)+2*Npad);
tempim(Npad+1:Npad+Nraw(1),Npad+1:Npad+Nraw(2)) = imraw;

tempim(1:Npad,Npad+1:Npad+Nraw(2)) = ones(Npad,1)*imraw(1,:);
tempim(Npad+1+Nraw(1):end,Npad+1:Npad+Nraw(2)) = ones(Npad,1)*imraw(end,:);

tempim(Npad+1:Npad+Nraw(1),1:Npad) = imraw(:,1)*ones(1,Npad);
tempim(Npad+1:Npad+Nraw(1),Npad+1+Nraw(2):end) = imraw(:,end)*ones(1,Npad);

imraw = tempim;

N = size(imraw);

% imshow (imraw/max(max(imraw)))
% return

fftim = fft2(imraw);

% 
% fftim = abs(double(real(fftim))) ;
% imshow(fftim/(0.01*max(max(fftim))))
% fftim(1:10,1:10)
% return

% 
% 
% fftim = fft2(imraw);

%create gaussian filter
% 1. fill quadrants of matrix with i^2 + j^2; ATTENTION: FFT stops at
% N-1, so you need two different quadrants and think about scaling with 2*N
% indexmatquadrant =(((0:((N(1)/2)-1) )/(N(1))).^2)'*ones(1,(N(2)/2)) + ones((N(1)/2),1)*(((0:((N(2)/2)-1))/(N(2))).^2);
indexmatquadrant_a = (((0:((N(1)/2)-1) )/(N(1))).^2)'*ones(1,(N(2)/2)) + ones((N(1)/2),1)*(((0:((N(2)/2)-1))/(N(2))).^2);
indexmatquadrant_b = fliplr((((0:((N(1)/2)-1) )/(N(1))).^2)'*ones(1,(N(2)/2)) + ones((N(1)/2),1)*(((1:((N(2)/2)))/(N(2))).^2));
indexmatquadrant_c = flipud((((1:((N(1)/2)) )/(N(1))).^2)'*ones(1,(N(2)/2)) + ones((N(1)/2),1)*(((0:((N(2)/2)-1))/(N(2))).^2));
indexmatquadrant_d = rot90(((( (1:((N(1)/2)))/(N(1))).^2)'*ones(1,(N(2)/2)) + ones((N(1)/2),1)*(( (1:((N(2)/2)))/(N(2))).^2)),2);


%2. fill quadrants of filtermatrix with gaussian, mean: 0, sigma:
% epsilon; ATTENTION: FFT in matlab is defined with factor 2*pi in the
% frequencies 

filterquadrant_a = exp(-( 2*pi^2*(epsin^2) )*indexmatquadrant_a);
filterquadrant_b = exp(-(2*pi^2*(epsin^2))*indexmatquadrant_b);
filterquadrant_c = exp(-( 2*pi^2*(epsin^2) )*indexmatquadrant_c);
filterquadrant_d = exp(-(2*pi^2*(epsin^2))*indexmatquadrant_d);



% 3. fill the whole filtermatrix
filter = [filterquadrant_a , filterquadrant_b ; filterquadrant_c , filterquadrant_d];

% imshow(filter/max(max(filter)))
% return

%4. inverse discrete Fourier transform of filtered image

imnew = ifft2(fftim.*filter,'symmetric');

% imshow(imnew/max(max(imnew)))
% return

% get original size; ATTENTION: has a small imaginary part due to rounding
% errors

%imnew = real(imnew(51:Nstart(1)+50,51:Nstart(2)+50));
imnew = real(imnew);
imsub = imraw - imnew;


imsub = imsub(Npad+1:Npad+Nraw(1),Npad+1:Npad+Nraw(2));
imnew = imnew(Npad+1:Npad+Nraw(1),Npad+1:Npad+Nraw(2));

