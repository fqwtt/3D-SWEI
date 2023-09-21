function transVec = getTransVector(dataPath)
RPath = 'I:\data\11162\tempR.txt';
tPath = 'I:\data\11162\tempt.txt';
seperateRt(dataPath, RPath, tPath);
R = load(RPath);
R = reshape(R', 3, 3, []);
R = permute(R, [2, 1, 3]);
t = load(tPath);

% Calculate World Coordinates
boardHeight = 11;
boardWidth = 4;
squareSize = 5; %mm
[Wx, Wy, Wz] = calWorldPoints(boardHeight, boardWidth, squareSize);

% Caluculate Camera Coordinates
length = size(R, 3);
CamCords = zeros(length, boardHeight * boardWidth, 3);
for i = 1 : length
    for j = 1 : boardHeight * boardWidth
        CamCords(i, j, :) = R(:,:,i) * [Wx(j); Wy(j); Wz(j)] + t(i, :)';
    end
end

xyz0 = zeros(boardHeight * boardWidth,3);
direction = zeros(boardHeight * boardWidth,3);
for j = 1 : boardHeight * boardWidth
    [xyz0(j,:), direction(j,:)] = fitLine3D(CamCords(:, j, 1), ...
        CamCords(:, j, 2), CamCords(:, j, 3), false);
end

transVec = mean(direction, 1);
transVec = transVec / norm(transVec);

end