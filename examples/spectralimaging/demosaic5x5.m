function cube = demosaic5x5(frame)

I=frame;
rowstart = 4;
colstart = 1;

D = single(zeros(216,409,16));


% This map translates band number to (row,col) position of the band
map = [1 1; 1 2; 1 3; 1 4; 1 5
       2 1; 2 2; 2 3; 2 4; 2 5
       3 1; 3 2; 3 3; 3 4; 3 5
       4 1; 4 2; 4 3; 4 4; 4 5;
       5 1; 5 2; 5 3; 5 4; 5 5]-1;

% Demosaic
for k = 1:25
    for row=1:(size(D,1)-1)
        for col=1:(size(D,2)-1)
            D(row,col,k) = I(rowstart+5*(row-1)+map(k,1),colstart+5*(col-1)+map(k,2));
        end
    end
    
    end
    cube= D;
end