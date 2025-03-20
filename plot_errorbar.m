function [h1,h2]=myerrorbar(x,y,std,drawline,linecolor,wid,markersize,flip)

if isvector(x) && ~isrow(x)
    x = x';
elseif ~isrow(x)
    error('x has to be a row')
end
if isvector(y) && ~isrow(y)
    y = y';
elseif ~isrow(y)
    error('y has to be a row')
end
if isvector(std) && ~isrow(std)
    std = std';
elseif ~isrow(std)
    error('std has to be a row')
end
if length(x) ~= length(y)
    error('x and y must be same length')
end
if ~exist('markersize','var') || isempty(markersize)
    markersize = 5;
end
if ~exist('wid','var') || isempty(wid)
    wid = 2;
end
if ~exist('drawline','var') || isempty(drawline)
    drawline = 1;
end
if ~exist('linecolor','var') || isempty(linecolor)
    linecolor = 'b';
end
if ~exist('flip','var') || isempty(flip)
    flip=0;
end

if flip==0
    h1 = plot([x;x],[(y-std);(y+std)]); hold on;
else
    h1 = plot([x-std;x+std],[(y);(y)]); hold on;
end
if drawline==1
    h2 = plot(x,y,'-o','MarkerSize',markersize);
elseif drawline==0
    h2 = plot(x,y,'o','MarkerSize',markersize);
end

if exist('linecolor','var')
    set(h1,'Color',linecolor)
    set(h2,'Color',linecolor,'MarkerFaceColor',linecolor)
end
if exist('wid','var')
    set(h1,'LineWidth',wid)
    set(h2,'LineWidth',wid)
end

h1 = h1(1);
h2 = h2(1);