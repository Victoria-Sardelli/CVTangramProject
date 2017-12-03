function findshapeTest

clear all;
close all;
clc;

squareM = [-0.5 -0.5 0.5 0.5; 0.5 -0.5 -0.5 0.5; 0 0 0 0];
triM = [-1/4 -1/4 3/4; -3/4 1/4 1/4; 0 0 0];
paralM = [-3*sqrt(2)/4 sqrt(2)/4 3*sqrt(2)/4 -sqrt(2)/4;...
           sqrt(2)/4   sqrt(2)/4 -sqrt(2)/4  -sqrt(2)/4;...
           0           0         0           0];
paralMAddPt = [0         0          sqrt(2)/2 -sqrt(2)/2;...
               sqrt(2)/4 -sqrt(2)/4 0         0;...
               0         0          0         0];
%paralM = [paralM paralMAddPt];
paralM(2,:) = -paralM(2,:);
       
squareRot = [0 89];
triRot = [0 359];
paralRot = [0 179];

%I = rgb2gray(im2double(imread('TangramSquare.png')));
%houghShapeFigure(I, [124 128], 1, [0 359], 45, triM, 'Triangle');
%Iedges = edge(I, 'Canny');
%figure(1), imshow(Iedges);

I = rgb2gray(im2double(imread('TestShapes.png')));
%Iparal = I(:,600:end);
%Isquare = padarray(I(:,1:220), [100 200], 0, 'both');
%houghShapeFigure(Isquare, 200, 5, 135, 45, 0, paralM, 'Paralellogram S:200 R:0 No flip');

%houghShapeFigure(I, 200, 5, 45, 45, 0, paralM, 'Paralellogram S:200 R:45 No flip');
%houghShapeFigure(I, 200, 5, 0, 45, 0, triM, 'Triangle S:200 R:0');
%houghShapeFigure(I, 200, 10, 0, 45, 0, squareM, 'Square S:200 R:0');

%houghShapeFigure(I, [100 210], 10, squareRot, 45, 0, squareM, 'Square');
%houghShapeFigure(I, [190 210], 5, triRot, 45, 0, triM, 'Triangle');
findShapeFigure(I, [190 210], 5, paralRot, 45, 1, paralM, 'Paralellogram');
%houghShapeFigure(I, [190 210], 5, paralRot, 45, 1, paralM, 'Paralellogram');

end

function gradImFigure(I)
    %h = fspecial('gaussian', 5, 2);
    h = ones(5) / 25;
    I = imfilter(I,h);
    [Igx, Igy] = imgradientxy(I, 'sobel');
    Ig = atand(Igy ./ Igx);
    [nr,nc] = size(Ig);
    figure;
    hold on
    colormap hsv
    colorbar
    pcolor([Ig nan(nr,1); nan(1,nc+1)]);
    shading flat;
    set(gca, 'ydir', 'reverse');
    hold off
end

function findShapeFigure(I, sizeRange, sizeInc, rotRange, rotInc, flip, shapeTP, fTitle)
    [origin, rotation, scale, flip, n] =  imfindshapes(I, sizeRange, 'SizeIncrement', sizeInc,...
        'RotationRange', rotRange, 'RotationIncrement', rotInc,...
        'Flip', flip, 'ShapeTemplate', shapeTP);
    %% Initiate figure
    f = figure;
    im = imshow(I);
    hold on;
    %im.AlphaData = 0.5;
    ax = gca;
    for i = [1:n]
        plotshape(shapeTP, origin(i,:), scale(i), rotation(i), 'Flip', flip(i), 'Axis', ax);
    end
end

function houghShapeFigure(I, sizeRange, sizeInc, rotRange, rotInc, flip, shapeTP, fTitle)
    %% Get hough shape transform
    [houghtrans, rotRange, sizeRange, flipRange] = houghshapes(I, sizeRange, 'SizeIncrement', sizeInc,...
        'RotationRange', rotRange, 'RotationIncrement', rotInc,...
        'Flip', flip, 'ShapeTemplate', shapeTP);
    
    %% Get and print maxes
    maxVal = max(max(max(max(max(houghtrans,[],4),[],3),[],5)))

    maxI = find(houghtrans == maxVal);
    [y,x,r,s,f] = ind2sub(size(houghtrans), maxI)

    valForI = ones(length(y), 1) * NaN;
    for i = [1:length(y)]
        valForI(i) = houghtrans(y(i), x(i), r(i), s(i), f(i));
        disp(['Val at (' num2str(y(i)) ', ' num2str(x(i)) ', ' num2str(rotRange(r(i))) ', ' ...
            num2str(sizeRange(s(i))) ', ' num2str(flipRange(f(i))) '): ' num2str(valForI(i))]);
    end
    
    %% Get dimension sizes for rotation and size
    rSize = size(houghtrans,3);
    sSize = size(houghtrans,4);
    
    %% Set/Calculate position values
    axPos = [0.13 0.31 0.615 0.615]; % Constant value for position
    axLblH = 0.05;
    uiUnits = 'normalized';
    compHRatio = [1 1 3/4 0]; % Slider Button Label Padding
    compBaseH = (axPos(2)-axLblH) / (2*compHRatio(1) + compHRatio(2) + ...
        4*compHRatio(3) + 5*compHRatio(4));
    sliderH = compBaseH * compHRatio(1);
    buttonH = compBaseH * compHRatio(2);
    lblH = compBaseH * compHRatio(3);
    paddingH = compBaseH * compHRatio(4);
    
    sliderW = axPos(3);
    sLblW = sliderW / sSize;
    rLblW = sliderW / rSize;
    buttonW = sliderW / 2;
        
    sliderX = axPos(1);
    buttonX = sliderX + (sliderW - buttonW)/2;
    
    buttonY = paddingH;
    rMLblY = buttonY + buttonH + paddingH;
    rSliderY = rMLblY + lblH;
    rLblY = rSliderY + sliderH;
    sMLblY = rLblY + lblH + paddingH;
    sSliderY = sMLblY + lblH;
    sLblY = sSliderY + sliderH;
    
    
    %% Initiate figure
    f = figure;
    im = imshow(I);
    hold on;
    im.AlphaData = 0.5;
    axI = gca;
    axI.Position = axPos;
    
    axHT = axes('position', axPos);
    htFrame = imagesc(houghtrans(:,:,1,1,1)); %, [0 size(shapeTP, 2)]
    htFrame.AlphaData = 0.5;
    colormap(axHT, hot);
    colorbar;
    title(['HST: ' fTitle]);
    axHT.Position = axPos;
    axHT.Color = 'none';
    axHT.PlotBoxAspectRatio = axI.PlotBoxAspectRatio;
    hold off;
    linkaxes([axI,axHT],'xy');
    
    % Clear varibles from setting up figure
    figVarList = {'im', 'axI', 'axHT'};
    clear(figVarList{:});
    
    % Set up controls
    bgcolor = [0.9400    0.9400    0.9400];
    
    % max = 0.31
    % rSliderY + (2*sliderH) + (3*lblH)
    % (2*sliderH) + (4*lblH) <= 0.31
    % lblH * (2*23/15 + 4) <= 0.31
    % lblH * (46/15 + 75/15)
    % lblH * (131/15) = lblH*8.7333
    % lblH <= 0.0355
    
    % lblH:sliderH:buttonH = 15:23:23
    % OR -- lblH:sliderH:buttonH = 15:20:20 = 3:4:4 -- EASIER
    % axY = 0.31
    % uiH <= axY
    % 4lblH/3 = sliderH
    % uiH = 2sliderH + 4lblH + buttonH + 5padY
    %     = 2sliderH + 4lblH + 5padY
    %     = lblH(19/4) + 5padY
    
%     uiUnits = 'normalized';
%     sliderX = 0.13;
%     sliderW = 0.615;
%     sliderH = 0.0544;
%     lblH = 0.0355;
%     
%     % VALUES BASED ON OTHER VALUES
%     rSliderY = lblH;
%     sSliderY = rSliderY + sliderH + (2*lblH);
%     rMLblY = rSliderY - lblH;
%     sMLblY = sSliderY - lblH;
%     sLblY = sSliderY + sliderH;
%     sLblW = sliderW / sSize;
%     rLblY = rSliderY + sliderH;
%     rLblW = sliderW / rSize;

    
    
    ri = 1;
    si = 1;
    fi = 1;
    
    % Scale Slider
    if sSize > 1
        sSlider = uicontrol('Parent',f,'Style','slider','units', uiUnits,...
            'Position',[sliderX,sSliderY,sliderW,sliderH],'value',si,...
            'min',1, 'max',sSize, 'SliderStep', [1/sSize 1/sSize]);
        for i = [0:sSize-1]
            lblX = sliderX + (i * sLblW);
            lblStr = num2str(min(sizeRange) + (i * sizeInc));
            uicontrol('Parent',f,'Style','text','units', uiUnits,...
                'Position',[lblX,sLblY,sLblW,lblH],...
                'String',lblStr,'BackgroundColor',bgcolor);
        end
        slbl = uicontrol('Parent',f,'Style','text','units', uiUnits,...
            'Position',[sliderX,sMLblY,sliderW,lblH],'String',...
            ['Scale: ' num2str((si-1)*sizeInc + min(sizeRange))],...
            'BackgroundColor',bgcolor);
        sSlider.Callback = @(es,ed) sCallback(es,ed);
    else
        slbl = uicontrol('Parent',f,'Style','text','units', uiUnits,...
            'Position',[sliderX,sSliderY,sliderW,lblH],'String',...
            ['Scale: ' num2str((si-1)*sizeInc + min(sizeRange))],...
            'BackgroundColor',bgcolor);
    end
    
    % Rotation Slider
    if rSize > 1
        rSlider = uicontrol('Parent',f,'Style','slider','units', uiUnits,...
            'Position',[sliderX,rSliderY,sliderW,sliderH],'value',ri,...
            'min',1, 'max', rSize, 'SliderStep', [1/rSize 1/rSize]);
        for i = [0:rSize-1]
            lblX = sliderX + (i * rLblW);
            lblStr = num2str(min(rotRange) + (i * rotInc));
            uicontrol('Parent',f,'Style','text','units', uiUnits,...
                'Position',[lblX,rLblY,rLblW,lblH],...
                'String',lblStr,'BackgroundColor',bgcolor);
        end
        rlbl = uicontrol('Parent',f,'Style','text','units', uiUnits,...
            'Position',[sliderX,rMLblY,sliderW,lblH],'String',...
            ['Rotation: ' num2str((ri-1)*rotInc + min(rotRange))],...
            'BackgroundColor',bgcolor);
        rSlider.Callback = @(es,ed) rCallback(es,ed);
    else
        rlbl = uicontrol('Parent',f,'Style','text','units', uiUnits,...
            'Position',[sliderX,rSliderY,sliderW,sliderH],'String',...
            ['Rotation: ' num2str((ri-1)*rotInc + min(rotRange))],...
            'BackgroundColor',bgcolor);
    end
    
    % Flip toggle button
    if flip
        fToggle = uicontrol('Parent',f,'Style','togglebutton',...
            'units', uiUnits,'Position',[buttonX buttonY buttonW buttonH],...
            'String','Flip','Min',1,'Max',2,'Value',fi,...
            'Callback', @(es,ed) fCallback(es,ed));
    end

    
    % Clear varibles from setting up ui
    uiVarList = {'bgcolor', 'uiUnits', 'sliderX', 'sliderW', 'sliderH',...
        'lblH', 'rSliderY', 'sSliderY', 'rMLblY', 'sMLblY', 'sLblY',...
        'sLblW', 'rLblY', 'rLblW', 'axPos'};
    clear(uiVarList{:});

    %% Callback functions
    function rCallback(es,ed)
        ri = round(es.Value);
        set(htFrame, 'CData', houghtrans(:,:,ri,si,fi));
        set(rlbl, 'String', ['Rotation: '...
            num2str((ri-1)*rotInc + min(rotRange))]);
    end
    function sCallback(es,ed)
        si = round(es.Value);
        set(htFrame, 'CData', houghtrans(:,:,ri,si,fi));
        set(slbl, 'String', ['Scale: '...
            num2str(min(sizeRange) + ((si-1) * sizeInc))]);
    end
    function fCallback(es,ed)
        fi = es.Value;
        set(htFrame, 'CData', houghtrans(:,:,ri,si,fi));
        if fi == 1
            fToggle.String = 'Flip';
        else
            fToggle.String = 'Unflip';
        end
    end
end

% function rCallback = getRCallback(imsc, lbl, hought, rStep)
%     function rCallbackFunc(source, event)
%         ri = round((source.Value + (rStep/2)) / rStep);
%         set(imsc, 'CData', hought(:,:,ri,1));
%         set(lbl, 'String', ['Rotation: ' num2str((ri-1)*rStep)]);
%     end
%     rCallback = @rCallbackFunc;
% end



% triR = rotz(90) * triM;
% %triR = triR + (5/3);
% triR(3, :) = [];
% polyM = (triR * 20) + repmat([50;50],[1 3]);
% polyM = reshape(polyM, [1 length(polyM(:))]);
% 
% I = zeros(101,101,3);
% Ia = insertShape(I, 'filledPolygon', polyM, 'color', 'white', 'Opacity', 1);
% imagesc(Ia);