function [Data] = importData()

path = 'E:\nju307_wt\SWI_m\1129\';
folderNum = 8;

%% 导入R数据

Data.R = cell(folderNum,1);
for i = 1:folderNum
    filepath = strcat(path,string(i),'\','fileRtBack.txt');
    tempFile = importdata(filepath);
    n = length(tempFile)/2;
    R_w2c = zeros(3,3,n);

    for j = 1:n
        R_w2c(:,:,j) = reshape(tempFile(2*j-1,1:9),3,3)';
    end
    Data.R{i} = R_w2c;
end

%% 导入mat数据

Data.Mat=cell(folderNum,1);
for i = 1:folderNum
    filepath = strcat(path,string(i),'\','test\');
    n = length(dir(strcat(filepath,'*.mat')));
    mat = zeros(200,44,3,n);
    for j = 1:n
        load(strcat(filepath,string(j),'.mat'))
        mat(:,:,:,j) = res1;
    end
    Data.Mat{i} = mat;
end

%% 导入t2w数据

tempFile = importdata('fileRt_t2w_L7_3.txt');
Data.t2w = reshape(tempFile(1,1:9),3,3)';

%% 导入Rw2Lw数据

tempFile = importdata('fileRt_R2L_hexa.txt');
Data.Rw2Lw = reshape(tempFile(1,1:9),3,3)';

%% 导入两个板子在两个相机下的空间位置关系

% Data.w2c = zeros(3,3,2,folderNum-1);
% for i = 1:folderNum-1
%     filepath = strcat(path,string(i),'\','Rt.txt');
%     tempFile = importdata(filepath);
%     Data.w2c(:,:,1,i) = reshape(tempFile(1,1:9),3,3)';
%     Data.w2c(:,:,2,i) = reshape(tempFile(3,1:9),3,3)';
% end

%% 导入所有相机与相机1的空间位置关系

Data.cs2c1=zeros(3,3,folderNum);
filepath = strcat('cs2c1.txt');
cs2c1 = importdata(filepath);
for i=1:folderNum
    Data.cs2c1(:,:,i)=reshape(cs2c1(i,1:9),3,3)';
end

%% 导入ks数据

Data.ks=cell(folderNum,1);
for i = 1:folderNum
    filepath = strcat(path,string(i),'\');
    load(strcat(filepath,'ks.mat'));
    Data.ks{i} = ks;
end

end


