function [movieParam,detectionParam] = makeParams( utracksettingsfile )

% movieParam
%              'imageDir',...
%              'filenameBase',...
%              'firstImageNum',...
%              'lastImageNum',...
%              'digits4Enum'

tld_ = regexp( utracksettingsfile,'(.*)(?<=\\)','match'); tld_ = tld_{1};

movieParam.('filenameBase') = 'fov1_';
movieParam.('imageDir') = strcat( tld_, 'ImageData\\' );
movieParam.('firstImageNum') = 1;
movieParam.('lastImageNum') = numel(dir( movieParam.('imageDir') ))-2;
movieParam.('digits4Enum') = 5;

% detectionParam
params = {'psfSigma',...
            'visual',...
            'doMMF',...
            'bitDepth',...
            'alphaLocMax',...
            'numSigmaIter',...
            'integWindow',...
            'maskLoc',...
            'calcMethod',...
            'testAlphaR',...
            'testAlphaA',...
            'testAlphaD',...
            'testAlphaF'};

fid = fopen( utracksettingsfile );
all_text = textscan( fid, '%s' );

for param = params
   detectionParam.(param{1}) = str2double(all_text{1}{find(cellfun( @(x) numel(x)>0, strfind(all_text{1},param{1}) )==1)+1});
end

% Workaround in case of linebreaks (folders should always use //)
detectionParam.('maskLoc') = join( all_text{1}(find(cellfun( @(x) numel(x)>0, strfind(all_text{1},'maskLoc') ) == 1)+1 : find(cellfun( @(x) numel(x)>0, strfind(all_text{1},'calcMethod') ) == 1)-1) );
detectionParam.('maskLoc') = detectionParam.('maskLoc'){1};
if eq( exist( detectionParam.('maskLoc') ), 0 ); detectionParam.('maskLoc') = strcat( tld_, 'mask.tif' ); end
if eq( exist( detectionParam.('maskLoc') ), 0 ); fprintf('No suitable mask file located!\n'); return; end;
detectionParam.('calcMethod') = all_text{1}{find(cellfun( @(x) numel(x)>0, strfind(all_text{1},'calcMethod') ) == 1)+1};

end