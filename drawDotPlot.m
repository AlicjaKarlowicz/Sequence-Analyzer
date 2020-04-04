function [image] = drawDotPlot(filMat)

% draw dot mat as image 
figure;
image = imagesc(filMat);
colormap(1-gray); 
title('Filtered dot matrix');
xlabel('sequence 2'); ylabel('sequence 1');
ax = gca;
set(ax,'XAxisLocation','top','Units','pixels');


end

