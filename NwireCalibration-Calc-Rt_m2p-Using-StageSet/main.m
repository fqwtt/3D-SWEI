clear;
clc;
close all;

%% Settings

dataPath = 'I:\data\11163\';
% The length of filey and filez can be changed (at least 1 of each)
filey = strings(1,1);
% !!!the first pair of R and t in the last of filey will be
% considered standard (phantom and axes coordination systems considered
% overlapped) 
% filey(1) = 'fileRt_w2c_y3.txt';
% filey(2) = 'fileRt_w2c_y2.txt';
filey(1) = 'fileRt_thetay.txt';
filez = strings(1,1);
% filez(1) = 'fileRt_w2c_z3.txt';
% filez(2) = 'fileRt_w2c_z2.txt';
filez(1) = 'fileRt_thetaz.txt';

% Store temporary data.  Can be anywhere.
RPath = 'I:\data\11163\tempR.txt';
tPath = 'I:\data\11163\tempt.txt';

% Whether use translation vector as the normal vector in circle fitting
useTransVectorY = false;
useTransVectorZ = false;
diry = [];
dirz = [];
if useTransVectorY == true
    pathy = [dataPath, 'fileRt_y.txt'];
    diry = getTransVector(pathy);
    if diry(2) < 0
        diry = -diry;
    end
end
if useTransVectorZ == true
    pathz = [dataPath, 'fileRt_z.txt'];
    dirz = getTransVector(pathz);
    if dirz(3) < 0
        dirz = -dirz;
    end
end

% if set to true -- axis Y is considered accurate, and axis Z is
% translated to meet axis Y; else axis Z will be considered accurate.
isYaccurate = true;
 
% true -- compute least squares of fitted circle centers to form axes lines
% false -- compute average of fitted circle centers and normals to form
% axes lines
LSQorAVRG = false;

% shift camera coordination system before calculate t_m2p, to reduce error.
shiftCamera = true;

%% Data Read and Process

% Calculate World Coordinates
boardHeight = 11;
boardWidth = 4;
squareSize = 6; %mm
[Wx, Wy, Wz] = calWorldPoints(boardHeight, boardWidth, squareSize);

ACenters_XY = [];
ACenters_XZ = [];
AN_XY = [];
AN_XZ = [];

figure;

for fileInd = 1 : length(filey)
    dataXZ = dataPath + filey(fileInd);
    seperateRt(dataXZ, RPath, tPath);
    R_XZ = load(RPath);
    R_XZ = reshape(R_XZ', 3, 3, []);
    R_XZ = permute(R_XZ, [2, 1, 3]);
    t_XZ = load(tPath);
    
    % Caluculate Camera Coordinates
    length_XZ = size(R_XZ, 3);
    CamCords_XZ = zeros(length_XZ, boardHeight * boardWidth, 3);
    for i = 1 : length_XZ
        for j = 1 : boardHeight * boardWidth
            CamCords_XZ(i, j, :) = R_XZ(:,:,i) * [Wx(j); Wy(j); Wz(j)] + t_XZ(i, :)';
        end
    end
    CamCords_XZ = permute(CamCords_XZ, [2, 1, 3]);
    
    % Fit Circles
    for i = 1 : boardHeight * boardWidth
        if useTransVectorY == true
            [Centers_XZ(i, :), N_XZ(i, :)] = fitCircle3D_norm(squeeze(CamCords_XZ(i, :, :)), diry);
        else
            [Centers_XZ(i, :), N_XZ(i, :)] = fitCircle3D(squeeze(CamCords_XZ(i, :, :)));
        end
    end
    ACenters_XZ = [ACenters_XZ; Centers_XZ];
    AN_XZ = [AN_XZ; N_XZ];
end

for fileInd = 1 : length(filez)
    dataXY = dataPath + filez(fileInd);
    seperateRt(dataXY, RPath, tPath);
    R_XY = load(RPath);
    R_XY = reshape(R_XY', 3, 3, []);
    R_XY = permute(R_XY, [2, 1, 3]);
    t_XY = load(tPath);
    
    % Caluculate Camera Coordinates
    length_XY = size(R_XY, 3);
    CamCords_XY = zeros(length_XY, boardHeight * boardWidth, 3);
    for i = 1 : length_XY
        for j = 1 : boardHeight * boardWidth
            CamCords_XY(i, j, :) = R_XY(:,:,i) * [Wx(j); Wy(j); Wz(j)] + t_XY(i, :)';
        end
    end
    CamCords_XY = permute(CamCords_XY, [2, 1, 3]);
    
    % Fit Circles
    for i = 1 : boardHeight * boardWidth
        if useTransVectorZ == true
            [Centers_XY(i, :), N_XY(i, :)] = fitCircle3D_norm( squeeze(CamCords_XY(i, :, :)),dirz );
        else
            [Centers_XY(i, :), N_XY(i, :)] = fitCircle3D( squeeze(CamCords_XY(i, :, :)) );
        end
    end
    ACenters_XY = [ACenters_XY; Centers_XY];
    AN_XY = [AN_XY; N_XY];
end

xlabel('x'); ylabel('y'); zlabel('z');

%% Find Rotation Axis

if LSQorAVRG == true
    [xyz0_XY, direction_XY] = fitLine3D(ACenters_XY(:,1), ACenters_XY(:,2),...
        ACenters_XY(:,3), false);
    [xyz0_XZ, direction_XZ] = fitLine3D(ACenters_XZ(:,1), ACenters_XZ(:,2),...
        ACenters_XZ(:,3), false);
else
    xyz0_XY = mean(ACenters_XY, 1);
    direction_XY = mean(AN_XY, 1);
    xyz0_XZ = mean(ACenters_XZ, 1);
    direction_XZ = mean(AN_XZ, 1);
end

if direction_XY(3) < 0
    direction_XY = -direction_XY; 
end
if direction_XZ(2) < 0
    direction_XZ = -direction_XZ; 
end

%% Get transform matrices as calibration results

% Convert camera coordinating to rotation-axis coordinating.
% The last parameter shifts the original point of rotation-axis coordinate 
% system along y direction, to make the original point lies on the center
% of the rotating platform's upper surface
[R_c2a, t_c2a] = getCam2Axes(xyz0_XZ, direction_XZ, xyz0_XY, direction_XY, 35, isYaccurate);

% Because t_c2a amplifies error in t_m2p, we shift camera coordinate system
% to change t_c2a to zeros
if shiftCamera == true
    for i = 1 : 3
        t_XY(:,i) = t_XY(:,i) + t_c2a(i);
        t_XZ(:,i) = t_XZ(:,i) + t_c2a(i);
    end
    t_c2a = [0; 0; 0];
end

% Convert marker coordinating to phantom coordinating.
% Marker to Cam * Cam to Phantom = Marker to Phantom
% ***axes and phantom coordination must be the same (both platforms set to 0)
% in the selected frame***
R_m2p = R_XZ(:, :, 1) * R_c2a ;
t_m2p = R_XZ(:, :, 1) * t_c2a + t_XZ(1, :)';
               
% Reproject marker points to camera coordination for testing purpose
P_test = zeros(boardHeight * boardWidth, 3);
for i = 1 : boardHeight * boardWidth
    P_test(i, :) = R_m2p * [Wx(i); Wy(i); Wz(i)] + t_m2p;
end

%% Save Results
R_m2p_save = reshape(R_m2p', 1, 9);
t_m2p_save = t_m2p';
writematrix(R_m2p_save, [dataPath, 'R_m2p.txt'], 'Delimiter', ' ');
writematrix(t_m2p_save, [dataPath, 't_m2p.txt'], 'Delimiter', ' ');
