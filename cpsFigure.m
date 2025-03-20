function h = cpsFigure(width,height,num)

Idx = exist('num','var');
if Idx
    h = figure(num);
else
    h = figure;
end

%set(h,'DefaultAxesLooseInset',[0,0,0,0]);
Position = get(h,'Position');
Position(3) = width*Position(3);
Position(4) = height*Position(4);
set(h,'Position', Position,'color','w');
set(h,'PaperPosition', [0 0 width height]*10,'color','w');