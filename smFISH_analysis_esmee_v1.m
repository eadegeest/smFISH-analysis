%% Original scripts:
% * completefit_differentiation_290118.m
% * dotfitting_parameters.m
% * fftfilter_nowrap_fast_pad.m
% * findcandidates_3d_manual.m
% * improcess_2015.m
% 
% << Dependencies >>
% Main file: completefit_differentiation_290118.m
% 1. Select files (figures, stacks)
% 2. Choose dot fitting parameters [dotfitting_parameters.m]
% 3. Process images [improcess_2015.m]
%     * Background subtraction [fftfilter_nowrap_fast_pad.m]
%     * Denoising [fftfilter_nowrap_fast_pad.m]
% 4. Find candidates (smFISH dots) and select the right threshold manually [findcandidates_3d_manual.m] 
%     * Project 3D image stack to a 2D image using maximum projection
%     * Calculate noise from maximum projection
%     * Thresholding > find candidates
%     * Measure properties of each candidate
%     * Save candidates on image
% 5. (Dependent script [assignpeakstocells_2015] not included: assign dots to cells using DAPI image)

% Here: all necessary parts of the scripts are written here in one script,
%so you can easily go through all the steps yourself
%% Write your own wrap around to fit multiple files etc.

basedirs = {...
'/Users/esmeeadegeest/Documents/Projects/TRICK/smFISH analysis in Matlab/Figures/Figure 1/'
};

%% Choose dot fitting parameters (the ones you need only)

    % thresholding parameter ( pxl > thresh*noise)
    thresh = 8;                 
    % dot fitting parameters
    sigma = 2.7; % PSF FWHM 
    %background filtering parameters
    filt = 2*sigma;
    denoise = 0;
    
%% Find figure(s)
basedir = basedirs{1};
mat = dir([basedir '*.tif']);
infile = [basedir mat(1).name(1:end-4)];

%% Process images: background subtraction and denoising
% << function [] = improcess_2015(infile,frame,bgfilt,denoise,plotswitch) >>
% "improcess_2015([infile '.tif'],[1 inf],filt,0,0);"
infile = [infile '.tif']; % file name of tif stack
frame = [1 inf]; % vector containg first and last slice to be processed e.g. [1 22]
bgfilt = filt; % all structures larger than this value will be removed (should be at least 3 times FWHM of dot)
denoise = 0; % all structures smaller than this value will be removed (should be zero for very dense signals)
plotswitch = 1; % set to 1 to see output of processing, otherwise set to 0

tifname = [infile(1:end-4) '_filt' num2str(bgfilt) '_denoise' num2str(denoise) '.tif'];

% load movie
info = imfinfo(infile);
if frame(2) > numel(info)
    frame(2) = numel(info);
end

framenum = frame(1):frame(2);
framenum(1);     

im = imread(infile,framenum(1),'Info',info);

%% PT1 - Background subtraction
% << function[imnew,imsub] = fftfilter_nowrap_fast_pad(imraw,epsin) >>
% "[imnew,imsub] = fftfilter_nowrap_fast_pad(im,bgfilt);"      
imraw = im; % raw image
epsin = bgfilt; % parameter epsilon, all structures on length scales > epsin will
                    % be removed (= spatial frequencies below 1/epsin)             
imraw = double(imraw);
Nraw = size(imraw);
Npad = 50; 
    tempim = zeros(size(imraw,1)+2*Npad,size(imraw,2)+2*Npad); % 1024 > 1124; 1124*1124
    tempim(Npad+1:Npad+Nraw(1),Npad+1:Npad+Nraw(2)) = imraw; % 51:1074,
    tempim(1:Npad,Npad+1:Npad+Nraw(2)) = ones(Npad,1)*imraw(1,:);
    tempim(Npad+1+Nraw(1):end,Npad+1:Npad+Nraw(2)) = ones(Npad,1)*imraw(end,:);
    tempim(Npad+1:Npad+Nraw(1),1:Npad) = imraw(:,1)*ones(1,Npad);
    tempim(Npad+1:Npad+Nraw(1),Npad+1+Nraw(2):end) = imraw(:,end)*ones(1,Npad);
    imraw = tempim;
    N = size(imraw);
    
    fftim = fft2(imraw);
    
% Create gaussian filter
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
filterquadrant_a = exp(-(2*pi^2*(epsin^2))*indexmatquadrant_a);
filterquadrant_b = exp(-(2*pi^2*(epsin^2))*indexmatquadrant_b);
filterquadrant_c = exp(-(2*pi^2*(epsin^2))*indexmatquadrant_c);
filterquadrant_d = exp(-(2*pi^2*(epsin^2))*indexmatquadrant_d);
    % 3. fill the whole filtermatrix
filter = [filterquadrant_a , filterquadrant_b ; filterquadrant_c , filterquadrant_d];
imshow(filter,[])
    % 4. inverse discrete Fourier transform of filtered image
imnew = ifft2(fftim.*filter,'symmetric');
imnew = real(imnew);
imsub = imraw - imnew;
imsub = imsub(Npad+1:Npad+Nraw(1),Npad+1:Npad+Nraw(2));
imnew = imnew(Npad+1:Npad+Nraw(1),Npad+1:Npad+Nraw(2));

imnew_b = imnew; % store background image under a new name 'imnew_b'
    % OUTPUT parameters
    % imnew: the background image (= only structures > epsin)
    % imsub: raw image - background image    

    
%% PT2 - Denoising
% << function[imnew,imsub] = fftfilter_nowrap_fast_pad(imraw,epsin) >> 
% "[imout,imsubb] = fftfilter_nowrap_fast_pad(imsub,denoise);"
imraw = imsub;
epsin = denoise; % parameter epsilon, all structures on length scales > epsin will
                    % be removed (= spatial frequencies below 1/epsin) 
             
imraw = double(imraw);
Nraw = size(imraw);
Npad = 50;
    
    tempim = zeros(size(imraw,1)+2*Npad,size(imraw,2)+2*Npad); % 1024 > 1124; 1124*1124
    tempim(Npad+1:Npad+Nraw(1),Npad+1:Npad+Nraw(2)) = imraw; % 51:1074,
    tempim(1:Npad,Npad+1:Npad+Nraw(2)) = ones(Npad,1)*imraw(1,:);
    tempim(Npad+1+Nraw(1):end,Npad+1:Npad+Nraw(2)) = ones(Npad,1)*imraw(end,:);
    tempim(Npad+1:Npad+Nraw(1),1:Npad) = imraw(:,1)*ones(1,Npad);
    tempim(Npad+1:Npad+Nraw(1),Npad+1+Nraw(2):end) = imraw(:,end)*ones(1,Npad);
    imraw = tempim;
    N = size(imraw);
    
    fftim = fft2(imraw);
    
% Create gaussian filter
    % 1. fill quadrants of matrix with i^2 + j^2; ATTENTION: FFT stops at
% N-1, so you need two different quadrants and think about scaling with 2*N
% indexmatquadrant =(((0:((N(1)/2)-1) )/(N(1))).^2)'*ones(1,(N(2)/2)) + ones((N(1)/2),1)*(((0:((N(2)/2)-1))/(N(2))).^2);
indexmatquadrant_a = (((0:((N(1)/2)-1) )/(N(1))).^2)'*ones(1,(N(2)/2)) + ones((N(1)/2),1)*(((0:((N(2)/2)-1))/(N(2))).^2);
indexmatquadrant_b = fliplr((((0:((N(1)/2)-1) )/(N(1))).^2)'*ones(1,(N(2)/2)) + ones((N(1)/2),1)*(((1:((N(2)/2)))/(N(2))).^2));
indexmatquadrant_c = flipud((((1:((N(1)/2)) )/(N(1))).^2)'*ones(1,(N(2)/2)) + ones((N(1)/2),1)*(((0:((N(2)/2)-1))/(N(2))).^2));
indexmatquadrant_d = rot90(((( (1:((N(1)/2)))/(N(1))).^2)'*ones(1,(N(2)/2)) + ones((N(1)/2),1)*(( (1:((N(2)/2)))/(N(2))).^2)),2);
   % 2. fill quadrants of filtermatrix with gaussian, mean: 0, sigma:
% epsilon; ATTENTION: FFT in matlab is defined with factor 2*pi in the
% frequencies 
filterquadrant_a = exp(-(2*pi^2*(epsin^2))*indexmatquadrant_a);
filterquadrant_b = exp(-(2*pi^2*(epsin^2))*indexmatquadrant_b);
filterquadrant_c = exp(-(2*pi^2*(epsin^2))*indexmatquadrant_c);
filterquadrant_d = exp(-(2*pi^2*(epsin^2))*indexmatquadrant_d);
    % 3. fill the whole filtermatrix
filter = [filterquadrant_a , filterquadrant_b ; filterquadrant_c , filterquadrant_d];
imshow(filter,[])
    % 4. inverse discrete Fourier transform of filtered image
imnew = ifft2(fftim.*filter,'symmetric');

imnew = real(imnew);
imsub = imraw - imnew;
imsub = imsub(Npad+1:Npad+Nraw(1),Npad+1:Npad+Nraw(2));
imnew = imnew(Npad+1:Npad+Nraw(1),Npad+1:Npad+Nraw(2));
    % OUTPUT parameters
imout = imnew;     % imout: background subtracted and denoised image
imsubb = imsub;    % imsubb: background subtracted image - denoised image          

%% Show eventual image and save it
    if plotswitch
        %plotim = imout/max(max(imout));
        %mean(plotim(plotim > 0))
        %imshow(plotim,[0 5*mean(plotim(plotim > 0))]);
        imshow(imout/max(max(imout)))
        pause(0.1)
    end
    imwrite(uint16(imout),tifname,'WriteMode','append') % save filtered image
    
%% Find candidates ('peaks')
% << function[peaks] = findcandidates_3d_manual(infile,frame,thresh,sigma,filtsigma) >>
% "peaks = findcandidates_3d_manual([infilefilt '.tif'],[1 inf],thresh,sigma,filt);" 
% find peak candidates in filtered image; thresholding with fixed threshold "thresh"
mat = dir([basedir '*_denoise0*.tif']);
infilefilt = [basedir mat(1).name(1:end-4)]; 
infile = [infilefilt '.tif'];
frame = [1 inf];
thresh = thresh;
sigma = sigma;
filtsigma = filt;

% round sigma to an integer
sigma = ceil(sigma);
% load movie
info = imfinfo(infile);
if frame(2) > numel(info)
    frame(2) = numel(info);
end
im = imread(infile,1,'Info',info);
%
bwmask = ones(size(im));
% get width and height
width = info(1).Width
height = info(1).Height
numpxl = width*height;
peaks = [];
maxim =zeros(height,width);
bwmat = zeros(frame(2)-frame(1)+1,height,width);
rawmat = zeros(frame(2)-frame(1)+1,height,width);

countframe = 0;

% Project 3D image stack to a 2D image: maximum projection
for framenum = frame(1):frame(2)
    framenum
    countframe = countframe + 1;

    % get raw image
    imraw = imread(infile,framenum,'Info',info);
    imraw = double(imraw);
    
    rawmat(countframe,:,:) = imraw;
    tempim = imraw;
    
    % get the maximum value of each pixel compared to all frames (in case of a stack)
    ind = find(tempim > maxim);
    if ~isempty(ind)
        maxim(ind) = tempim(ind);
    end
end

% Calculate noise from maximum projection
noisepow = (abs(fft(maxim(1:end))).^2)/numpxl; % power of DFT = power of noise
noise = sqrt(mean(noisepow(round(numpxl/(2*filtsigma)):round(numpxl/2)))); % sq.rt. of average power of noise for the 'frequencies' higher than 
% noise in higher frequency range, noise with spatial properties smaller than FWHM 

% Thresholding + measure properties of each candidate
done = 0;
tempthresh = thresh % thresholding parameter ( pxl > thresh*noise)

while done < 1

imthresh = max([tempthresh*noise,1]); % max. intensity value: tempthresh*noise or 1?
    
bwmat = logical(rawmat > imthresh); % which pixels have a value larger than threshold? #Candidate particles

    CC = bwconncomp(bwmat,6); % Find connected components in binary image. desired connectivity: three-dimensional six-connected neighborhood
    N = CC.NumObjects; % number of connected components = # of candidates

    bwprops = regionprops(CC,{'Centroid','Area','PixelIdxList','BoundingBox'}); % Measure properties of image regions (black-white)
    grayprops = regionprops(CC,rawmat,{'WeightedCentroid','MinIntensity','MeanIntensity','MaxIntensity'}); % Measure properties of image regions (gray)
  
    peaks = [cat(1,bwprops.Centroid),cat(1,grayprops.WeightedCentroid)]; % concatenate info
    
     % bring x,y,z in order
    tempvec = peaks(:,2);
    peaks(:,2) = peaks(:,3);
    peaks(:,3) = tempvec;

    tempvec = peaks(:,5);
    peaks(:,5) = peaks(:,6);
    peaks(:,6) = tempvec;
    
    grayim = maxim/quantile(maxim(1:end),0.9999);
    outim = zeros(size(grayim,1),size(grayim,2),3);
    outim(:,:,1) = grayim;
    outim(:,:,2) = grayim;
    outim(:,:,3) = grayim;
    imshow(outim)
    hold on
    plot(peaks(:,5),peaks(:,4),'or','MarkerSize',3)
    hold off
    
    waitforbuttonpress
    
    key = get(gcf,'CurrentCharacter');
    
    if key == 'n'
        done = 1
    elseif key == ','
        tempthresh = max([1,tempthresh-1])
    elseif key == '.'
        tempthresh = tempthresh+1
    
    end
end

peaks = [];
count = 0;

for i = 1:N % for each candidate
%i=1;
    bbox = bwprops(i).BoundingBox;
    if bbox(5) > 0
        count = count + 1;
        if mod(count,1000)==0
            count
        end
bbox(1:3) = round(bbox(1:3));

        indxmin = max([1,bbox(2)-1]);
        indxmax = min([size(bwmat,1),indxmin+bbox(5)+2]);

        indymin = max([1,bbox(1)-1]);
        indymax = min([size(bwmat,2),indymin+bbox(4)+2]);

        indzmin = max([1,bbox(3)-1]);
        indzmax = min([size(bwmat,3),indzmin+bbox(6)+2]);

        smallrawmat = rawmat([indxmin:indxmax],[indymin:indymax],[indzmin:indzmax]);
        smallbwmat = bwmat([indxmin:indxmax],[indymin:indymax],[indzmin:indzmax]);
        
        indpxl = bwprops(i).PixelIdxList;
        [indx,indy,indz] = ind2sub(size(bwmat),indpxl);
        
        indx = indx - indxmin + 1;
        indy = indy - indymin + 1;
        indz = indz - indzmin + 1;

        subim = zeros(size(smallbwmat));
        
        for k = 1:size(indx)
            subim(indx(k),indy(k),indz(k)) = 1;        
        end

        subim = imdilate(subim,ones(2,2,2));

        %get perimeter of that region
        bw2 = bwperim(subim,26);

        grayperim = bw2.*smallrawmat;
        Nperim = sum(bw2(:));

        bg1 = sum(grayperim(:))/Nperim;
        
        %exclude other peaks from background calculation
        bw2 = bw2.*(ones(size(smallbwmat))-smallbwmat);

        grayperim = bw2.*smallrawmat;
        Nperim = sum(bw2(:));

        bg2 = sum(grayperim(:))/Nperim;

        peaks(count,:) = [bwprops(i).Centroid,grayprops(i).WeightedCentroid,bwprops(i).Area,...
            grayprops(i).MinIntensity,grayprops(i).MeanIntensity,grayprops(i).MaxIntensity,bbox(4),bbox(6),bbox(5),bg1,bg2];
    end
end

% bring x,y,z in order
tempvec = peaks(:,2);
peaks(:,2) = peaks(:,3);
peaks(:,3) = tempvec;

tempvec = peaks(:,5);
peaks(:,5) = peaks(:,6);
peaks(:,6) = tempvec;


% save check image
grayim = maxim/quantile(maxim(1:end),0.9999);
outim = zeros(size(grayim,1),size(grayim,2),3);
outim(:,:,1) = grayim;
outim(:,:,2) = grayim;
outim(:,:,3) = grayim;

for i = 1:size(peaks,1)
    
    if peaks(i,end) > 0
        outim(round(peaks(i,4)),round(peaks(i,5)),1) = 1;
        outim(round(peaks(i,4)),round(peaks(i,5)),2) = 0;
        outim(round(peaks(i,4)),round(peaks(i,5)),3) = 0;
    else
        outim(round(peaks(i,4)),round(peaks(i,5)),1) = 1;
        outim(round(peaks(i,4)),round(peaks(i,5)),2) = 0;
        outim(round(peaks(i,4)),round(peaks(i,5)),3) = 0;
    end
    
end

imwrite(double(outim),[infile(1:end-4) '_check.tif'])

size(peaks)

save( [infilefilt '.pkc'],'peaks','-ASCII');  
