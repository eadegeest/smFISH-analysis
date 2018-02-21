% write your own wrap around to fit multiple files etc.

basedirs = {...
'/Users/esmeeadegeest/Documents/Projects/TRICK/smFISH analysis in Matlab/Figures/Figure 1/HighContrast/Fig_1/'
};

% dot fitting parameters
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
%

for i = 1:size(basedirs,1)

    basedir = basedirs{i};
    
    % background filtering
    mat = dir([basedir '*.tif']);
    

%     matdapi = dir([basedir '/dapi/*.tif'])

%     
    
    for k = 1:size(mat,1)

        infile = [basedir mat(k).name(1:end-4)];

        % high pass fourier filter
        improcess_2015([infile '.tif'],[1 inf],filt,0,0);

    end
%     

    % peak finding
    mat = dir([basedir '*_denoise0*.tif']);

    for k = 1:size(mat,1)

            infilefilt = [basedir mat(k).name(1:end-4)];           
   
            % find peak candidates
 
            % thresholding with fixed threshold "thresh"
            peaks = findcandidates_3d_manual([infilefilt '.tif'],[1 inf],thresh,sigma,filt);                       
           %   peaks = findcandidates_3d_manual_withrefine([infilefilt '.tif'],[1 inf],thresh,sigma,filt,areathresh,delflood);                       
%            peaks = findcandidates_3d_automatic([infilefilt '.tif'],[1 inf],thresh,sigma,filt,areathresh,delflood);                       
            save( [infilefilt '.pkc'],'peaks','-ASCII');

    end
    
    % %manual segmentation

%    matdapi = dir([basedirs{i} 'dapi/*.tif'])
%     matdapi = dir([basedirs{i} '*/a594/*.tif'])

    %     for k = 1:size(matdapi,1)
%         infiledapi = [basedir '\dapi\' matdapi(k).name(1:end-4)]

%     for k = 1:size(matdapi,1)
% 
%             infiledapi = [basedirs{i} 'dapi/' matdapi(k).name(1:end-4)]           
% %             infile = [basedirs{i} mat(k).name(1:end-4)]
% 
%             %segment cells
%             
%              [X] = cellsegment_findcells_maxproject_interactive([infiledapi '.tif'],[0 inf],0.1);
% 
%             imwrite(uint8(X),[infiledapi '_Xmat.tif'],'tif');
%             close all
% 
%     end      


end

return





for i = 1:size(basedirs,1)

    basedir = basedirs{i};
    
    
%    assign dots to cells
    


    mata = dir([basedir '*c2*.pkc'])
    matb = [];
    matc = [];

    matXmat = dir([basedir '/dapi/*_Xmat.tif']);


    for k = 1:max([size(mata,1),size(matb,1),size(matc,1)])

            if ~isempty(matXmat)
                X = imread([basedir '/dapi/' matXmat(k).name]);
            else
                X = ones(1024,1024);
            end

            if ~isempty(mata)
                infilea = [basedir mata(k).name];   
                peaks = load(infilea);
                [outpkza] = assignpeakstocells_2015(peaks,nucradius,X);
                save( [infilea(1:end-4) '.pka'],'outpkza','-ASCII');
            end

            if ~isempty(matb)
                infileb = [basedir matb(k).name];  
                peaks = load(infileb);
                [outpkza] = assignpeakstocells_2015(peaks,nucradius,X);
                save( [infileb(1:end-4) '.pka'],'outpkza','-ASCII');
            end

            if ~isempty(matc)
                infilec = [basedir matc(k).name];  
                peaks = load(infilec);
                [outpkza] = assignpeakstocells_2015(peaks,nucradius,X);
                save( [infilec(1:end-4) '.pka'],'outpkza','-ASCII');
            end


    end

    % %collect data
    matXmat = dir([basedir '/dapi/*_Xmat.tif']);
    matgfpmat = dir([basedir '/tdtomato/*c1*.tif']);
    mata = dir([basedir '*c2*.pka'])
%     matb = dir([basedir '*a594*.pka'])
%     matc = dir([basedir '*tmr*.pka'])

    res = [];
    
    

    for k = 1:max([size(mata,1),size(matb,1),size(matc,1)])
        
        k                 

        if ~isempty(matXmat)
                X = imread([basedir '/dapi/' matXmat(k).name]);
        else
                X = ones(1024,1024);
        end

        if ~isempty(mata)
            infilea = [basedir mata(k).name];   
            peaksa = load(infilea);
        end

        if ~isempty(matb)
            infileb = [basedir matb(k).name];   
            peaksb = load(infileb);
        end

        if ~isempty(matc)
            infilec = [basedir matc(k).name];   
            peaksc = load(infilec);
        end

        
        [area,mintens] = getgfp([basedir '/tdtomato/' matgfpmat(k).name],X);
        
        
        
        M = max(max(X));

        for m = 1:M
            subim = (X == m);
            S = regionprops(subim,'Centroid');
                         res = [res; double(k),double(m), double(round(S.Centroid)), area(m), mintens(m), double(sum(peaksa(:,end)==m))];
%             res = [res; double(k),double(m), double(round(S.Centroid)), area(m), mintens(m), double(sum(peaksa(:,end)==m)), double(sum(peaksb(:,end)==m)),double(sum(peaksc(:,end)==m))];
        end

    %     inda = logical(peaksa(:,end) > -1);

    %     figure
    %     imshow(double(X)/double(max(max(X))))
    %     hold on
    %     plot(peaksa(inda,2),peaksa(inda,1),'or')
    %     plot(matb(:,1),matb(:,2),'om')
    %     plot(matc(:,1),matc(:,2),'og')

    end
    
    res
%      save([basedir '122115_leonardo_WT_2i_os.txt'],'res','-ASCII')

end


