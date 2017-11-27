function houghtrans =  imfindshapes(varargin)
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

    %% Convert image to grayscale
    I = rgb2gray(I);

    %% Get image edges
    Iedges = edge(I,'Canny');
    Iedges = im2double(Iedges);
    
    %% Get size and rotation ranges
    if (length(sizeRange) > 1)
        sizeRange = minSize:sizeInc:max(sizeRange);
    end
    if (length(rotRange) > 1)
        rotRange = min(rotRange):rotInc:max(rotRange);
    end
    
    %% Get shape template size
    shapeTPSize = max(sum(shapeTP));
    
    %% Create voting space and vote
    houghtrans = zeros(size(Iedges,1), size(Iedges,2),...
        length(sizeRange), length(rotRange));
    for si = [1:length(sizeRange)]
        s = sizeRange(si);
        resizedShape = imresize(shapeTP, s/shapeTPSize);
        for ri = [1:length(rotRange)]
            r = rotRange(ri);
            rotShape = imrotate(resizedShape, r);
            shapeEdges = im2double(edge(rotShape, 'Canny'));
            votingSpace = imfilter(Iedges, shapeEdges);
            houghtrans(:, :, si, ri) = votingSpace;
        end
    end

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
    parser.addParamValue('ShapeTemplate',[1 1; 1 1]); % @checkShapeTemplate
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
                'nonsparse','nonempty','positive','finite','vector'});
        else
            validateattributes(rotRange,{'numeric'},{'integer','nonnan', ...
                'nonsparse','nonempty','positive','finite','vector'});
        end
        if (length(rotRange) > 2)
            error(message('tangram:imfindshapes:unrecognizedSizeRange'));
        elseif (length(rotRange) == 2)
            if (rotRange(1) > rotRange(2))
                error(message('images:imfindshapes:invalidSizeRange'));
            elseif rotRange(2) > 359
                error(message('images:imfindshapes:invalidSizeRange'));
            end
        end
        
        tf = true;
    end

    function tf = checkFlip(flip)
        validateattributes(flip,{'logical'},{'nonnan', ...
            'nonsparse','nonempty','positive','finite','vector'},mfilename,'RADIUS_RANGE',2);
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
        
end