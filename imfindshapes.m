function [centers, rotations, scales, flip, n] =  imfindshapes(varargin)
% IMFINDSHAPES Uses hough transform to find shapes in an image
% houghtrans = imfindshapes(I,s) Returns voting space for squares of
% size s rotated from 0 to 90 degrees.
%
% See also HOUGH and HOUGHSHAPES.

    %% Get hough shape transform and ranges
    [houghtrans, rotRange, sizeRange, flipRange] =  houghshapes(varargin{:});
    
    %% Get max value and indices
    maxVal = max(max(max(max(max(houghtrans,[],4),[],3),[],5)));
    maxI = find(houghtrans == maxVal);

    %% Get subscripts of max value
    szHoughTrans = size(houghtrans);
    [y,x,r,s,f] = ind2sub(szHoughTrans, maxI);
    
    %% Get values to return
    centers = [x y];
    rotations = rotRange(r);
    scales = sizeRange(s);
    flip = flipRange(f);
    n = length(y);
end