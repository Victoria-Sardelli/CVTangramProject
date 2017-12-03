function [houghtrans, rotRange, sizeRange, flipRange] =  imfindshapes(varargin)
% IMFINDSHAPES Uses hough transform to find shapes in an image
% houghtrans = imfindsquares(I,s) Returns voting space for squares of
% size s rotated from 0 to 90 degrees.
%
% See also HOUGH.

%% Get inputs
parsedInputs = parseInputs(varargin{:});

I = parsedInputs.Image;
sizeRange = parsedInputs.SizeRange;
rotRange = parsedInputs.RotationRange;
sizeInc = parsedInputs.SizeIncrement;
rotInc = parsedInputs.RotationIncrement;
flip = parsedInputs.Flip;
shapeTP = parsedInputs.ShapeTemplate;
edgeThresh = parsedInputs.EdgeThreshold;

    %% Rotate shape template 180 degrees round origin
    %% ACTUALLY UNEEDED
    % Assuming the origin of the shape defined in the given points is
    % (0,0), inverting the points will find the value of
    % (x_c - x_{ij}, y_c - y_{ij}) that will be used when shifting points
    shapeM = -1 * shapeTP;
    shapeTP = shapeM;
    
    %% Get gradient values/direction at shape points
    % Plus 18 and mod 360 to keep values positive and in range [0,359]
    shapeGDir = mod(atan2d(shapeTP(1,:), shapeTP(2,:)) + 180, 360);
    
    %% Get max width of shape for padding and angles of edges...?
    % acosd
    pointCt = size(shapeM,2);
    maxLineL = 0;
    lineAngles = zeros(1, factorial(pointCt)/(2 * factorial(pointCt-2)));
    i = 1;
    for a = [1:pointCt-1]
        for b = [a+1:pointCt]
            lineL = sqrt(sum((shapeM(:,a) - shapeM(:,b)).^2));
            lineA = acosd(abs(shapeM(2,a) - shapeM(2,b))/lineL);
            if maxLineL < lineL
                maxLineL = lineL;
            end
            lineAngles(i) = lineA;
            i = i + 1;
        end
    end
    imPadding = ceil(maxLineL * max(sizeRange));
    lineAngles = unique(lineAngles);
    

    %% Convert image to grayscale, if needed, and get image size
    if (size(I,3) > 1)
        I = rgb2gray(I);
    end
    [Iheight, Iwidth] = size(I);
    
    %% Get image gradient, and pad with zeros
    [Igx, Igy] = imgradientxy(I, 'sobel');
    %Ig = mod(atan2d(Igx,Igy) + 180, 360);
    Ig = atand(Igy./Igx);
    IgPadded = padarray(Ig, [imPadding imPadding], NaN, 'post');

    %% Get image edges, and pad with zeros
    Iedges = edge(I,'sobel');
    Iedges = im2double(Iedges);
    Iedges = imdilate(Iedges, ones(3)); % Dilate to thicken edges
    %Iedges = imrotate(Iedges,180);
    IedgesPadded = padarray(Iedges, [imPadding imPadding], 0, 'post');
    
    %% Get size and rotation ranges
    if (length(sizeRange) > 1)
        sizeRange = min(sizeRange):sizeInc:max(sizeRange);
    end
    if (length(rotRange) == 2)
        rotRange = min(rotRange):rotInc:max(rotRange);
        %possRot = unique(bsxfun(@minus,unique(Ig),shapeG))
    end
    if flip
        flipRange = [ 0 1 ];
    else
        flipRange = [ 0 ];
    end
    
    %% Create voting space and vote
    %houghtrans = zeros(size(Iedges,1), size(Iedges,2),...
    %    length(sizeRange), length(rotRange));
    houghtrans = zeros(Iheight, Iwidth,...
        length(rotRange), length(sizeRange), length(flipRange));
    for si = [1:length(sizeRange)]
        s = sizeRange(si);
        for ri = [1:length(rotRange)]
            r = rotRange(ri);
            for f = flipRange
                % Create shape transformation matrix
                tM = rotz(r) * diag([s*(-1)^f s s]);
                % Transform shape
                tShapeM = tM * shapeTP;
                % Get gradient directions for transformed shape
                %tShapeGDir = mod(abs((360*f) - shapeGDir) + r, 360);
                tShapeGDir = atand(tShapeM(2,:)./tShapeM(1,:)); % For now
                %tShapeGDir = mod(shapeGDir + r, 360);
                
                %spSz = ceil(sqrt(size(tShapeM,2)));
                %figure;
                
                vsG = zeros(size(IgPadded));
                vsE = zeros(size(IedgesPadded));
                for pi = [1:size(tShapeM,2)] % Loop over points
                    p = round(tShapeM(1:2,pi)); % Probably a more accurate option than rounding p
                    p = flipud(p); % Have to flip x and y for circle shift
                    %p(1) = -p(1); % Invert y because of image axis
                    pGDir = tShapeGDir(pi);
                    vsG = vsG + circshift(mod((IgPadded + 90) - (pGDir + 90),180) <= 1, p);
                    %votingSpace = circshift(IgPadded == pGDir, p);
%                     vsEtemp = circshift(IedgesPadded, p);
%                     h = subplot(spSz, spSz, pi);
%                     imshow(cat(3, Iedges, zeros(size(I)), zeros(size(I))));
%                     hold on;
%                     i = imagesc(vsEtemp(1:Iheight, 1:Iwidth));
%                     i.AlphaData = 0.5;
%                     plotshape(shapeTP,[Iwidth/2, Iheight/2],s,r,'Axis',h);
%                    vsE = vsE +vsEtemp;
                    vsE = vsE + circshift(IedgesPadded, p);
                end
%                 figure;
%                 imshow(cat(3, Iedges, zeros(size(I)), zeros(size(I))));
%                 hold on;
                i = imagesc(vsE(1:Iheight, 1:Iwidth));
                i.AlphaData = 0.5;
                    
                votingSpace = vsG + vsE;
                votingSpace = votingSpace(1:Iheight, 1:Iwidth);
                houghtrans(:, :, ri, si, f+1) = votingSpace;
            end
            
%             tShapeM(3,:) = []; % Get rid of z-coordinates
%             for p = tShapeM % Loop over points
%                 IedgesShift = circshift(IedgesPadded, round(p));
%                 % There is probably a more accurate option than rounding p
%                 IedgesShift = IedgesShift(1:Iheight, 1:Iwidth);
%                 %f = figure, subplot(2,1,1), imshow(Iedges),
%                 %subplot(2,1,2), imshow(IedgesShift);
%                 %close(f);
%                 houghtrans(:, :, ri, si) = houghtrans(:, :, ri, si) + IedgesShift;
%             end

            % smooth voting space
            %h = fspecial('gaussian',floor(s / 2),2);
            %houghtrans(:, :, ri, si) = imfilter(houghtrans(:, :, ri, si), h);
        end
    end
    
    % NEED TO FACTOR IN FILP AND GRADIENTS
    % imgradientxy(I, 'sobel');

end



% Function to parse inputes, copied from infindcircles
function parsedInputs = parseInputs(varargin)

narginchk(2,Inf);

persistent parser;

if (isempty(parser))
    parser = inputParser();

    parser.addRequired('Image',@checkImage);
    parser.addRequired('SizeRange',@checkSizeRange);
    %parser.addRequired('RotRange',@checkRotRange);
    %parser.addParamValue('Method','phasecode',@checkMethod);
    %parser.addParamValue('ObjectPolarity','bright');
    parser.addParamValue('RotationRange',[0 90],@checkRotRange);
    parser.addParamValue('Flip',0,@checkFlip);
    %parser.addParamValue('ShapeTemplate',[1 1; 1 1]); % @checkShapeTemplate
    parser.addParamValue('ShapeTemplate',[-0.5 -0.5 0.5 0.5;...
        -0.5 0.5 0.5 -0.5; 0 0 0 0]); % @checkShapeTemplate
    parser.addParamValue('SizeIncrement',1);
    parser.addParamValue('RotationIncrement',1);
    parser.addParamValue('EdgeThreshold',[]);
end

% Parse input, replacing partial name matches with the canonical form.
if (nargin > 3) % If any name-value pairs are given
  varargin(4:end) = images.internal.remapPartialParamNames({'RotationRange',...
      'SizeIncrement', 'RotationIncrement',...
      'Flip', 'ShapeTemplate', 'EdgeThreshold'}, varargin{4:end});
end

parser.parse(varargin{:});
parsedInputs = parser.Results;

validateSizeRange(); % If Smin and Smax are the same, setS = Smin
validateRotRange(); % If Rmin and Rmax are the same then set R = Rmin.
validateShapeTP();
    
    function tf = checkImage(A)
        allowedImageTypes = {'uint8', 'uint16', 'double', 'logical', 'single', 'int16'};
        validateattributes(A,allowedImageTypes,{'nonempty',...
            'nonsparse','real'},mfilename,'A',1);
        N = ndims(A);
        if (isvector(A) || N > 3)
            error(message('tangram:imfindshapes:invalidInputImage'));
        elseif (N == 3)
            if (size(A,3) ~= 3)
                error(message('tangram:imfindshapes:invalidImageFormat'));
            end
        end
        tf = true;
    end

    function tf = checkSizeRange(sizeRange)
        if (isscalar(sizeRange))
            validateattributes(sizeRange,{'numeric'},{'nonnan', ...
                'nonsparse','nonempty','positive','finite','vector'});
        else
            validateattributes(sizeRange,{'numeric'},{'integer','nonnan', ...
                'nonsparse','nonempty','positive','finite','vector'});
        end
        if (length(sizeRange) > 2)
            error(message('tangram:imfindshapes:unrecognizedSizeRange'));
        elseif (length(sizeRange) == 2)
            if (sizeRange(1) > sizeRange(2))
                error(message('images:imfindshapes:invalidSizeRange'));
            end
        end
        
        tf = true;
    end

    function tf = checkRotRange(rotRange)
        if (isscalar(rotRange))
            validateattributes(rotRange,{'numeric'},{'nonnan', ...
                'nonsparse','nonempty','finite','vector'});
        else
            validateattributes(rotRange,{'numeric'},{'integer','nonnan', ...
                'nonsparse','nonempty','finite','vector'});
        end
        %if (length(rotRange) > 2)
        %    error(message('tangram:imfindshapes:unrecognizedSizeRange'));
        %else
        if max(rotRange) > 359
                error(message('images:imfindshapes:invalidSizeRange'));
        end
        if (length(rotRange) == 2)
            if (rotRange(1) > rotRange(2))
                error(message('images:imfindshapes:invalidSizeRange'));
            end
        end
        
        tf = true;
    end

    function tf = checkFlip(flip)
        validateattributes(flip,{'logical', 'numeric'}, {'<=', 1});
        tf = true;
    end

    function validateSizeRange
        if (length(parsedInputs.SizeRange) == 2)
            if (parsedInputs.SizeRange(1) == parsedInputs.SizeRange(2))
                parsedInputs.SizeRange = parsedInputs.SizeRange(1);
            end
        end
    end

    function validateRotRange
        if (length(parsedInputs.RotationRange) == 2)
            if (parsedInputs.RotationRange(1) == parsedInputs.RotationRange(2))
                parsedInputs.RotationRange = parsedInputs.RotationRange(1);
            end
        end
    end

    function validateShapeTP
        if (size(parsedInputs.RotationRange,1) == 2)
            parsedInputs.RotationRange(3,:) = zeros(1,...
                size(parsedInputs.RotationRange,2));
        end
    end
        
end