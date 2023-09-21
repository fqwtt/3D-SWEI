function [px, py, pz] = calWorldPoints(boardHeight, boardWidth, squareSize)
px = zeros(boardHeight * boardWidth, 1);
py = zeros(boardHeight * boardWidth, 1);
pz = zeros(boardHeight * boardWidth, 1);
count = 0;
for i = 1: boardHeight
    for j = 1 : boardWidth
        count = count + 1;
        px(count) = (2 * (j-1) + mod(i-1, 2)) * squareSize;
        py(count) = (i-1) * squareSize;    
    end
end
end