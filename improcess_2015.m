function [] = improcess_2015(infile,frame,bgfilt,denoise,plotswitch)
% infile : file name of tif stack
% frame  : vector containg first and last slice to be processed e.g. [1 22]
% bgfilt : all structures larger than this value will be removed (should be at least 3 times FWHM of dot)
% denoise : all structures smaller than this value will be removed (should be zero for very dense signals)
% plotswitch : set to 1 to see output of processing, otherwise set to 0

tifname = [infile(1:end-4) '_filt' num2str(bgfilt) '_denoise' num2str(denoise) '.tif'];


% load movie
info = imfinfo(infile);
if frame(2) > numel(info)
    frame(2) = numel(info);
end




for framenum = frame(1):frame(2)
    framenum
            
    im = imread(infile,framenum,'Info',info);

    %first bg subtraction, then denoising    
    [imnew,imsub] = fftfilter_nowrap_fast_pad(im,bgfilt);        
    [imout,imsubb] = fftfilter_nowrap_fast_pad(imsub,denoise);
     
    if plotswitch
        %plotim = imout/max(max(imout));
        %mean(plotim(plotim > 0))
        %imshow(plotim,[0 5*mean(plotim(plotim > 0))]);
        imshow(imout/max(max(imout)))
        pause(0.1)
    end
    
    
    imwrite(uint16(imout),tifname,'WriteMode','append')


end

