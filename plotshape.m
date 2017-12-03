function plotshape(varargin)
    % Parse inputs
    parsedInputs = parseInputs(varargin{:});
    shapepts = parsedInputs.Shape;
    or = parsedInputs.Origin;
    sc = parsedInputs.Scale;
    rot = parsedInputs.Rotation;
    f = parsedInputs.Flip;
    ax = parsedInputs.Axis;
    color = parsedInputs.Color;
    lineW = parsedInputs.LineWidth;
    lineSty = parsedInputs.LineStyle;

    % Get X & Y points for shape
    scM = diag([sc*(-1)^f sc sc]);
    rotM = rotz(rot);
    shapepts = rotM * scM * shapepts;
    shapepts(1,:) = shapepts(1,:) + or(1);
    shapepts(2,:) = shapepts(2,:) + or(2);
    X = [shapepts(1,:) shapepts(1,1)];
    Y = [shapepts(2,:) shapepts(2,1)];
    
    % Plot shape
    plot(ax,X,Y,'Color',color,'LineWidth',lineW,'LineStyle',lineSty);
end

% Function to parse inputes, copied from infindcircles
function parsedInputs = parseInputs(varargin)

narginchk(2,Inf);

persistent parser;

if (isempty(parser))
    parser = inputParser();

    parser.addRequired('Shape',@checkShape);
    parser.addRequired('Origin',@checkOrigin);
    parser.addRequired('Scale',@checkScale);
    parser.addRequired('Rotation',@checkRotation);
    parser.addParamValue('Flip',0,@checkFlip);
    parser.addParamValue('Axis',gca);
    parser.addParamValue('Color','r');
    parser.addParamValue('LineWidth',0.5);
    parser.addParamValue('LineStyle','-');
end

% Parse input, replacing partial name matches with the canonical form.
if (nargin > 5) % If any name-value pairs are given
  varargin(5:end) = images.internal.remapPartialParamNames({'Flip',...
      'Axis', 'Color', 'LineWidth','LineStyle'}, varargin{5:end});
end

parser.parse(varargin{:});
parsedInputs = parser.Results;

    function tf = checkShape(shapeM)
        validateattributes(shapeM, {'numeric'}, {'nonempty', '2d',...
            'nrows', 3});
        tf = true;
    end
    function tf = checkOrigin(origin)
        validateattributes(origin, {'numeric'}, {'nonempty',...
            'numel', 2});
        tf = true;
    end
    function tf = checkScale(scale)
        validateattributes(scale, {'numeric'}, {'nonempty', 'nonzero',...
            'nonnan','numel', 1});
        tf = true;
    end
    function tf = checkRotation(rot)
        validateattributes(rot, {'numeric'}, {'nonempty','numel', 1,...
            'nonnegative', '<', 360});
        tf = true;
    end
    function tf = checkFlip(flip)
        validateattributes(flip, {'numeric', 'logical'}, {'nonempty',...
            'binary', 'numel', 1});
        tf = true;
    end
end