function mylabel(xlab, ylab, fontsize)

try   
labX = xlabel(xlab);
catch
end
try
labY = ylabel(ylab);
catch
end

if exist('fontsize','var')
    set(labX, 'FontSize', fontsize);
    set(labY, 'FontSize', fontsize);
end

end