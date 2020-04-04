function [filMat] = filterDotMat(dotMat,windowsize,errorlimit)

rownum = size(dotMat,1);
colnum = size(dotMat,2);

if windowsize == 0
    windowsize =1; 
end

% allocate space
filMat = zeros(rownum,colnum);

condition = windowsize - errorlimit; % min value of diagonal

if errorlimit > windowsize
    error('Error limit bigger or equal to window');
end

for i = 1:rownum - (windowsize-1)
    for j = 1:colnum - (windowsize-1)
        
        windowMat = dotMat(i:i+(windowsize-1),j:j+(windowsize-1));
        
        if sum(diag(windowMat)) >= condition % sum of diagonal of matrix
            filWin = diag(windowMat);
            filWin = diag(filWin);
            filMat(i:i+(windowsize-1),j:j+(windowsize-1)) = filWin;
        end
        
    end
    
end


end

