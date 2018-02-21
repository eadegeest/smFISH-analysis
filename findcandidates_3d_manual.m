function[peaks] = findcandidates_3d_manual(infile,frame,thresh,sigma,filtsigma)

% round sigma to an integer
sigma = ceil(sigma);

% load movie
info = imfinfo(infile);
if frame(2) > numel(info)
    frame(2) = numel(info);
end
   

im = imread(infile,1,'Info',info);
%imhandle = imshow(double(im)/max(max(double(im))));
%recthandle = imrect;
%bwmask = createMask(recthandle,imhandle);
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
% go through the designated frames

for framenum = frame(1):frame(2)
    framenum
    countframe = countframe + 1;

    % get raw image
    imraw = imread(infile,framenum,'Info',info);
    imraw = double(imraw);
    
    
    rawmat(countframe,:,:) = imraw;

    tempim = imraw;
    
    ind = find(tempim > maxim);
    if ~isempty(ind)
        maxim(ind) = tempim(ind);
    end
end
    
 

% calculate noise from max. projection
noisepow = (abs(fft(maxim(1:end))).^2)/numpxl;

noise = sqrt(mean(noisepow(round(numpxl/(2*filtsigma)):round(numpxl/2))));


done = 0;

tempthresh = thresh;

while done < 1

    imthresh = max([tempthresh*noise,1]);
    
    bwmat = logical(rawmat > imthresh);

    CC = bwconncomp(bwmat,6);
    N = CC.NumObjects


    bwprops = regionprops(CC,{'Centroid','Area','PixelIdxList','BoundingBox'});
    grayprops = regionprops(CC,rawmat,{'WeightedCentroid','MinIntensity','MeanIntensity','MaxIntensity'});

  
    peaks = [cat(1,bwprops.Centroid),cat(1,grayprops.WeightedCentroid)];
    


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

%     for i = 1:size(peaks,1)
% 
%         if peaks(i,end) > 0
%             outim(round(peaks(i,4)),round(peaks(i,5)),1) = 0;
%             outim(round(peaks(i,4)),round(peaks(i,5)),2) = 0;
%             outim(round(peaks(i,4)),round(peaks(i,5)),3) = 1;
%         else
%             outim(round(peaks(i,4)),round(peaks(i,5)),1) = 1;
%             outim(round(peaks(i,4)),round(peaks(i,5)),2) = 0;
%             outim(round(peaks(i,4)),round(peaks(i,5)),3) = 0;
%         end
% 
%     end

    imshow(outim)
    hold on
    plot(peaks(:,5),peaks(:,4),'or','MarkerSize',3)
    hold off
    
    waitforbuttonpress
    
    key = get(gcf,'CurrentCharacter');
    
    if key == 'n'
        done = 1;
    elseif key == ','
        tempthresh = max([1,tempthresh-1]);
    elseif key == '.'
        tempthresh = tempthresh+1;
    
    end
  
    %peaks
    %centroid x,y,z
    %weighted centroid x,y,z
    % area
    % min,mean,max pixel intensity
    % extension in x,y,z



end

 
% 


peaks = [];
count = 0;

for i = 1:N

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
% figure
% imshow(grayim)
% hold on
% plot(peaks(:,3),peaks(:,1),'xr')

