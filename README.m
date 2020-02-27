% February 19-25 %
%
% In MATLAB:
%
% For now, this file will function as a log of what we do during our in-person meetings

[movieParam,detectionParam] = makeParams( 'C:\u-track\#03Ch1\uTrackSettings.txt' )

% Edit: 2/26/2020
%
% You need to edit the uTrackSettings.txt file to set maskLoc to the location of mask.tif

% Except replace the C:\u-track\ with the location where you placed the data we downloaded today

% You can then execute the following, but I'd recommend placing a red dot at line 908 in detectSubResFeatures2D_StandAlone.m

detectSubResFeatures2D_StandAlone(movieParam,detectionParam,0,0)

% Once you're in there, you'll see all your workspace variables on the
% right. You can obtain an image of the current frame going through the
% loop with the following code:

figure; imagesc( imageRaw )

% See if you can use patch or plot to draw an ROI around the current object
% (featuresInfo contains the x- and y- coordinates and the uncertainty)

% Take a look at how the least squares, non-linear fitting is done. 
% You can find this on lines 320-329 of the similarily named file: detectSubResFeatures2D_V2

edit detectSubResFeatures2D_V2

% This routine tries to optimize an objective function, the difference
% between a theoretical image and the image of the point in the frame. You
% can find the objective function here:

edit fitNGaussians2DVarSigma

% See if you can output some quality-control parameters from the function,
% add any other constraints, distance transformation info...
% 
% We can talk about the rest next week!

% February 26-Mar 3
%
% We went through detectSubResFeatures2D_StandAlone (focusing on lines 315-673)

% Mainly, we looked at several of the important images that are generated in the pipeline:

figure; imagesc(bgMeanInteg) % Generated at line 426 using spatialMovAveBG.m
figure; imagesc(fImg) % Generated at line 453 using locmax2d.m
figure; imagesc(imageIntegF) % Generated at line 417 using filterGauss2D.m

% We looked at the "first filter" for passing potential molecules along to the next stage:

pValue = 1 - normcdf(localMaxAmp,bgMeanMaxF,bgStdMaxF); % At line 469

% if absBG
% bgAbsMeanMaxF = bgAbsMeanInteg(localMax1DIndx);
% bgAbsStdMaxF = bgAbsStdInteg(localMax1DIndx);
% pValueAbs = 1 - normcdf(localMaxAmp,bgAbsMeanMaxF,bgAbsStdMaxF);
% end
% %retain only those maxima with significant amplitude
% if absBG
% keepMax = find((pValue < alphaLocMax(iWindow)) & ...
% (pValueAbs < alphaLocMaxAbs));
% else
% keepMax = find(pValue < alphaLocMax(iWindow));
% end
% keepMax
% size(keepMax)
% pValue(keepMax)
% alphaLocMax
% figure; imagesc(imageIntegF)
% keepMax
% numel(keepMax)
% %construct cands structure
% if numLocalMax == 0 %if there are no local maxima
% cands = [];
% %                 emptyFrames = [emptyFrames; iImage+integWindow]; %#ok<AGROW>
% else %if there are local maxima
% %define background mean and status
% cands = repmat(struct('status',1,'IBkg',[],...
% 'Lmax',[],'amp',[],'pValue',[]),numLocalMax,1);
% %store maxima positions, amplitudes and p-values
% for iMax = 1 : numLocalMax
% cands(iMax).IBkg = bgMeanMax(iMax);
% cands(iMax).Lmax = [localMaxPosX(iMax) localMaxPosY(iMax)];
% cands(iMax).amp = localMaxAmp(iMax);
% cands(iMax).pValue = pValue(iMax);
% end
% end
% if absBG
% bgAbsMeanMaxF = bgAbsMeanInteg(localMax1DIndx);
% bgAbsStdMaxF = bgAbsStdInteg(localMax1DIndx);
% pValueAbs = 1 - normcdf(localMaxAmp,bgAbsMeanMaxF,bgAbsStdMaxF);
% end
% %retain only those maxima with significant amplitude
% if absBG
% keepMax = find((pValue < alphaLocMax(iWindow)) & ...
% (pValueAbs < alphaLocMaxAbs));
% else
% keepMax = find(pValue < alphaLocMax(iWindow));
% end
% localMaxPosX = localMaxPosX(keepMax);
% localMaxPosY = localMaxPosY(keepMax);
% localMaxAmp = localMaxAmp(keepMax);
% bgMeanMax = bgMeanMax(keepMax);
% pValue = pValue(keepMax);
% numLocalMax = length(keepMax);
% %construct cands structure
% if numLocalMax == 0 %if there are no local maxima
% cands = [];
% %                 emptyFrames = [emptyFrames; iImage+integWindow]; %#ok<AGROW>
% else %if there are local maxima
% %define background mean and status
% cands = repmat(struct('status',1,'IBkg',[],...
% 'Lmax',[],'amp',[],'pValue',[]),numLocalMax,1);
% %store maxima positions, amplitudes and p-values
% for iMax = 1 : numLocalMax
% cands(iMax).IBkg = bgMeanMax(iMax);
% cands(iMax).Lmax = [localMaxPosX(iMax) localMaxPosY(iMax)];
% cands(iMax).amp = localMaxAmp(iMax);
% cands(iMax).pValue = pValue(iMax);
% end
% end
% cands
% cands(2)
% figure; imagesc(imageIntegF)
% hold on; structfun( @(x) plot(x.Lmax,'kx'), cands )
% hold on; arrayfun( @(x) plot(x.Lmax,'kx'), cands )
% hold on; arrayfun( @(x) plot([x.Lmax(1),x.Lmax(2)],'kx'), cands )
% figure; imagesc(imageIntegF); hold on; arrayfun( @(x) plot([x.Lmax(1),x.Lmax(2)],'kx'), cands )
% cands
% cands.Lmax
% figure; plot(0,0); hold on; arrayfun( @(x) plot([x.Lmax(1),x.Lmax(2)],'kx'), cands )
% figure; plot(0,0); hold on; arrayfun( @(x) plot(x.Lmax(1),x.Lmax(2),'kx'), cands )
% figure; imagesc(imageIntegF); hold on; arrayfun( @(x) plot(x.Lmax(2),x.Lmax(1),'kx'), cands )
% figure; imagesc(imageIntegF); hold on; arrayfun( @(x) plot(cands(x).Lmax(2),cands(x).Lmax(1),'kx'), [1:numel(cands)] )
% arrayfun( @(x) text(cands(x).Lmax(2),cands(x).Lmax(1),x), [1:numel(cands)] )
% arrayfun( @(x) text(cands(x).Lmax(2),cands(x).Lmax(1),sprintf('%i',x) ), [1:numel(cands)] )
% cands(70)
% cands(27)
% %add the cands of the current image to the rest - this is done
% %for the raw images, not the integrated ones
% localMaxima(iImage+integWindow(iWindow)).cands = ...
% [localMaxima(iImage+integWindow(iWindow)).cands; cands];
% localMaxima
% localMaxima(1)
% localMaxima(2)
% iImage+integWindow(iWindow)
% localMaxima(3)
% integWindow
% localMaxima(iImage+integWindow(iWindow)).cands
% numImagesInteg
% dbquit
% detectSubResFeatures2D_StandAlone(movieParam,detectionParam,0,0)
% localMaxima
% localMaxima(4).cands
% localMaxima(6).cands
% localMaxima(9).cands
% localMaxima
% arrayfun( @(x) size( x.cands,1 ), localMaxima )
% figure; plot( arrayfun( @(x) size( x.cands,1 ), localMaxima ) )
% iImage
% iImage=3900;
% %store raw images in array
% imageRaw = NaN(imageSizeX,imageSizeY,1+2*integWindow(iWindow));
% for jImage = 1 : 1 + 2*integWindow(iWindow)
% if hasImageDir && imageExists(jImage+iImage-1)
% imageRaw(:,:,jImage) = double(imread([imageDir filenameBase ...
% enumString(imageIndx(jImage+iImage-1),:) '.tif']));
% elseif ~hasImageDir
% imageRaw(:,:,jImage) = double(channel.loadImage(imageIndx(jImage+iImage-1)));
% end
% end
% %apply mask
% imageRaw = imageRaw .* repmat(maskImage,[1 1 1+2*integWindow(iWindow)]);
% %replace zeros with NaNs
% %zeros result from cropping that leads to curved boundaries
% imageRaw(imageRaw==0) = NaN;
% %normalize images
% imageRaw = imageRaw / (2^bitDepth-1);
% %integrate images
% imageInteg = nanmean(imageRaw,3);
% %filter integrated image
% %         imageIntegF = filterGauss2D(imageInteg,min(1,psfSigma));
% imageIntegF = filterGauss2D(imageInteg,psfSigma);
% %use robustMean to get mean and std of background intensities
% %in this method, the intensities of actual features will look like
% %outliers, so we are effectively getting the mean and std of the background
% %account for possible spatial heterogeneity by taking a spatial moving
% %average
% %get integrated image background noise statistics
% [bgMeanInteg,bgStdInteg] = ...
% spatialMovAveBG(imageInteg,imageSizeX,imageSizeY);
% %calculate absolute background mean and std if supplied
% if absBG
% bgRaw = NaN(bgSizeX,bgSizeY,1+2*integWindow(iWindow));
% for jImage = 1 : 1 + 2*integWindow(iWindow)
% if imageExists(jImage+iImage-1)
% bgRaw(:,:,jImage) = double(imread([bgImageDir bgImageBase ...
% enumString(imageIndx(jImage+iImage-1),:) '.tif']));
% end
% end
% bgRaw(bgRaw==0) = NaN;
% bgRaw = bgRaw / (2^bitDepth-1);
% bgInteg = nanmean(bgRaw,3);
% bgAbsMeanInteg = nanmean(bgInteg(:))*ones(imageSizeX,imageSizeY);
% bgAbsStdInteg = nanstd(bgInteg(:))*ones(imageSizeX,imageSizeY);
% end
% background = [];
% %clear some variables
% clear imageInteg
% try
% %call locmax2d to get local maxima in filtered image
% fImg = locmax2d(imageIntegF,[3 3],1);
% %get positions and amplitudes of local maxima
% [localMaxPosX,localMaxPosY,localMaxAmp] = find(fImg);
% localMax1DIndx = find(fImg(:));
% %get background values corresponding to local maxima
% bgMeanInteg1 = bgMeanInteg;
% bgMeanMaxF = bgMeanInteg1(localMax1DIndx);
% bgStdInteg1 = bgStdInteg;
% bgStdMaxF = bgStdInteg1(localMax1DIndx);
% bgMeanMax = bgMeanRaw(localMax1DIndx);
% %calculate the p-value corresponding to the local maxima's amplitudes
% %assume that background intensity in integrated image is normally
% %distributed with mean bgMeanMaxF and standard deviation bgStdMaxF
% pValue = 1 - normcdf(localMaxAmp,bgMeanMaxF,bgStdMaxF);
% if absBG
% bgAbsMeanMaxF = bgAbsMeanInteg(localMax1DIndx);
% bgAbsStdMaxF = bgAbsStdInteg(localMax1DIndx);
% pValueAbs = 1 - normcdf(localMaxAmp,bgAbsMeanMaxF,bgAbsStdMaxF);
% end
% %retain only those maxima with significant amplitude
% if absBG
% keepMax = find((pValue < alphaLocMax(iWindow)) & ...
% (pValueAbs < alphaLocMaxAbs));
% else
% keepMax = find(pValue < alphaLocMax(iWindow));
% end
% localMaxPosX = localMaxPosX(keepMax);
% localMaxPosY = localMaxPosY(keepMax);
% localMaxAmp = localMaxAmp(keepMax);
% bgMeanMax = bgMeanMax(keepMax);
% pValue = pValue(keepMax);
% numLocalMax = length(keepMax);
% %construct cands structure
% if numLocalMax == 0 %if there are no local maxima
% cands = [];
% %                 emptyFrames = [emptyFrames; iImage+integWindow]; %#ok<AGROW>
% else %if there are local maxima
% %define background mean and status
% cands = repmat(struct('status',1,'IBkg',[],...
% 'Lmax',[],'amp',[],'pValue',[]),numLocalMax,1);
% %store maxima positions, amplitudes and p-values
% for iMax = 1 : numLocalMax
% cands(iMax).IBkg = bgMeanMax(iMax);
% cands(iMax).Lmax = [localMaxPosX(iMax) localMaxPosY(iMax)];
% cands(iMax).amp = localMaxAmp(iMax);
% cands(iMax).pValue = pValue(iMax);
% end
% end
% end
% figure; imagesc(imageIntegF);
% localMaxima
% localMaxima.cands(3999)
% localMaxima(3999).cands
% localMaxima(3995).cands
% figure; plot( arrayfun( @(x) size( x.cands,1 ), localMaxima ) )
% figure; plot( arrayfun( @(x) mean( x.cands.pValue ), localMaxima ) )
% figure; plot( arrayfun( @(x) mean( arrayfun(@(y) y.pValue, x.cands ) ), localMaxima ) )
% figure; plot( arrayfun( @(x) min( arrayfun(@(y) y.pValue, x.cands ) ), localMaxima ) )
% figure; plot( arrayfun( @(x) min( arrayfun(@(y) y.pValue, x.cands ), 'UniformOutput', false ), localMaxima ) )
% figure; plot( arrayfun( @(x) min( arrayfun(@(y) y.pValue, x.cands , 'UniformOutput', false ), localMaxima ) ) )
% figure; plot( arrayfun( @(x) min( arrayfun(@(y) y.pValue, x.cands ), [], 1 ), localMaxima ) )
% figure; plot( arrayfun( @(x) median( arrayfun(@(y) y.pValue, x.cands ) ), localMaxima ) )
% figure; plot( arrayfun( @(x) min( arrayfun(@(y) y.pValue, x.cands ) ), localMaxima ) )
% arrayfun(@(y) y.pValue, x.cands )
% arrayfun(@(y) y.pValue, localMaxima.cands )
% localMaxima(3600).cands
% localMaxima(3600).cands.pValue
% arrayfun( @(x) x, localMaxima(3600).cands.pValue )
% class( localMaxima(3600).cands.pValue )
% type( localMaxima(3600).cands.pValue )
% arrayfun(@(x) x.pValue, localMaxima(3600).cands )
% min(arrayfun(@(x) x.pValue, localMaxima(3600).cands ))
% min(arrayfun(@(x) x.pValue, localMaxima(end).cands ))
% min(arrayfun(@(x) x.pValue, localMaxima(3994).cands ))
% for i = 1:3994; min_pvalue(i) = min(arrayfun(@(x) x.pValue, localMaxima(i).cands )); end
% for i = 1:3994; min_pvalue(i) = min(min(arrayfun(@(x) x.pValue, localMaxima(i).cands ))); end
% min(arrayfun(@(x) x.pValue, localMaxima(i).cands ))
% i
% for i = 3:3994; min_pvalue(i) = min(min(arrayfun(@(x) x.pValue, localMaxima(i).cands ))); end
% figure; plot( min_pvalue )
% localMaxima(3800)
% localMaxima(3800).cands
% localMaxima(3800).cands.pValue
% figure; imagesc(imageIntegF);
% myrandmovie = rand(100,100,1000);
% size( myrandmovie )
% figure; for i = 1:1000; imagesc( myrandmovie(:,:,i) ); pause(0.03); end
% figure; for i = 1:1000; imagesc( myrandmovie(:,:,i) ); plot( rand(10,1), rand(10,1), 'ko'); pause(0.03); end
% figure; for i = 1:1000; imagesc( myrandmovie(:,:,i) ); hold on; plot( rand(10,1), rand(10,1), 'ko'); pause(0.03); hold off; end
% figure; for i = 1:1000; imagesc( myrandmovie(:,:,i) ); hold on; plot( 100*rand(10,1), 100*rand(10,1), 'ko'); pause(0.03); hold off; end
% figure; for i = 1:1000; imagesc( myrandmovie(:,:,i) ); hold on; plot( 100*rand(10,1), 100*rand(10,1), 'wo'); pause(0.03); hold off; end
% figure; for i = 1:1000; imagesc( myrandmovie(:,:,i) ); hold on; plot( 100*rand(10,1), 100*rand(10,1), 'rx'); pause(0.03); hold off; end
% figure; for i = 1:1000; imshow( myrandmovie(:,:,i) ); hold on; plot( 100*rand(10,1), 100*rand(10,1), 'rx'); pause(0.03); hold off; end
% figure; for i = 1:1000; imshow( myrandmovie(:,:,i), 'colormap', jet(200) ); hold on; plot( 100*rand(10,1), 100*rand(10,1), 'rx'); pause(0.03); hold off; end
% figure; for i = 1:1000; imshow( myrandmovie(:,:,i), 'colormap', parula(200) ); hold on; plot( 100*rand(10,1), 100*rand(10,1), 'rx'); pause(0.03); hold off; end
% figure; for i = 1:1000; imshow( myrandmovie(:,:,i), 'colormap', parula(200), 'initialmagnification', 400 ); hold on; plot( 100*rand(10,1), 100*rand(10,1), 'rx'); pause(0.03); hold off; end
% figure; for i = 1:1000; imshow( myrandmovie(:,:,i), [0,.1], 'colormap', parula(200), 'initialmagnification', 400 ); hold on; plot( 100*rand(10,1), 100*rand(10,1), 'rx'); pause(0.03); hold off; end
% figure; imagesc(imageIntegF);
% figure; imshow(imageIntegF,[0,0.001]);
% figure; imshow(imageIntegF,[0,0.001],'colormap', parula(200) );
