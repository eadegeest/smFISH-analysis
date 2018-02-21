%% Detect red blood cells using image segmentation
% Source: https://nl.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html

% Note: Use the background image created in 'PT1 - Background subtraction' of the
% script smFISH_analysis_esmee_v1.m. Make sure the image is smooth enough
% by using, for example, filt = 2*sigma.

% 1. Detect entire cell using binary gradient mask (find lines of high contrast)
% Use Zero-Cross ('zerocross') or Laplacian of Gaussian ('log') method for
% edge detection. Use fudgeFactor = 0.5 (default)
[~, threshold] = edge(imnew_b, 'log');
fudgeFactor = 0.5;
BWs = edge(imnew_b,'log', threshold * fudgeFactor);
close all
figure(1)
imshow(im,[]),title('original image')
figure(2)
imshow(BWs), title('binary gradient mask');

% 2. Dilate image (widen lines to remove linear gaps)
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
BWsdil = imdilate(BWs, [se90 se0]);
figure(3)
imshow(BWsdil), title('dilated gradient mask');

% 3. Fill interior gaps
BWdfill = imfill(BWsdil, 'holes');
figure(4), imshow(BWdfill);
title('binary image with filled holes');

% 4. Smoothen the segmented objects by eroding image twice with diamond 
% structure element.
seD = strel('diamond',1);
BWfinal = imerode(BWdfill,seD);
BWfinal = imerode(BWfinal,seD);
figure(5), imshow(BWfinal), title('segmented image');

% (5. Outline the cells)
BWoutline = bwperim(BWfinal);
Segout = im; 
Segout(BWoutline) = 255; 
figure(6), imshow(Segout,[]), title('outlined original image');