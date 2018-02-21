
  % cell segmentation freehand yes/no
    freehand = 1;

    % thresholding parameter ( pxl > thresh*noise)
    thresh = 8;                 

    areathresh = 7;
    delflood = 1;

    % dot fitting parameters
    sigma = 2.7; % PSF FWHM
    
    tolsigma = 0.1;
    gatesigma = 2.3;

    % nuclear size
    nucsize = 80;
    
    %background filtering parameters
    filt = 1*sigma;
    filtdapi = nucsize;
    
    denoise = 0;


    accur = 0.01; % accuracy for dot fitting (in pxl)
    maxiter = 100; % max iterations if no minimum is found

    %tracing parameters
    radius = 1;
    gap = 1;

    % tracing parameters
    mintracelength = 2; %tracelength > this value
    maxtracelength = inf; % tracelength < this value

    pxl = 0.130; % physical pxl size
    spacing = 0.2; % physical distance between layers



    % radius for assignment to cell
    nucradius  = 20;



    

    
    
    
    
    