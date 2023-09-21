%% 用颜色代表剪切波大小

%% 导入数据

% Data=importData();
XX=-10:0.1:10;
[XX,YY]=meshgrid(XX);

foldnum = 8; % 文件个数
cm=parula;   % color map

X=cell(foldnum,1);
Y=cell(foldnum,1);
Z=cell(foldnum,1);
points=cell(foldnum,1);
value=cell(foldnum,1);
pointss=cell(foldnum,1);
for i=1:foldnum
    n=length(Data.R{i});
    points{i}=zeros(n,3);
    
    %% 求trans的x轴在本相机下的向量
    for j=1:n
        temp=Data.R{i}(:,:,j)*Data.t2w;
        points{i}(j,:)=temp(:,1)';
    end
    
    %% 拟合平面并求出投影坐标
    [points{i}, b]=fitplane(points{i});
    
    %% 将点转移到xy平面

    t = [   b(1),        b(2),     b(1);
            b(2),       -b(1),     b(2);
          sum(b.^2),      0,       -1 ];
    t=t./sqrt(sum(t.^2));
    temp=t\points{i}';
    
  
    %% 求剪切波速度对应的颜色索引
    ks=min(abs(Data.ks{i}(:,1)),abs(Data.ks{i}(:,2)));
    x=temp(1,:)'.*ks;
    y=temp(2,:)'.*ks;
    x=[x;-x];
    y=[y;-y];
    x=x(x~=0);
    y=y(y~=0);
    W = fitellipse(x,y);
    % 绘制结果
    figure
    plot(x,y,'*'), axis equal
    hold on
    [x,y] = drawellipse(W);
    axis equal
    grid on
    value{i}=(ks/10).^index;
    value{i}=round(value{i}*255)+1;
    pointss{i}=zeros(length(x),3);
    pointss{i}(:,1)=x;
    pointss{i}(:,2)=y;
    pointss{i}=t*pointss{i}';

        %% 求trans的x轴在公共相机下的向量
    temp=Data.cs2c1(:,:,i)*pointss{i};
    r=(sum(temp.^2).^0.5); 
    pointss{i}=(temp);
end

close all
figure
for i=[1,2,6]
% for i=1:foldnum
    x=points{i}(:,1);
    y=points{i}(:,2);
    z=points{i}(:,3);
%     for j=1:length(x)
%         
%         plot3([x(j),-x(j)],[y(j),-y(j)],[z(j),-z(j)],'*','Color',cm(value{i}(j),:));
%         hold on;
%     end
    plot3(pointss{i}(1,:),pointss{i}(2,:),pointss{i}(3,:),'*');
    hold on;
end

radius=0.99;
[xx, yy, zz]=sphere(50); 
h2=surf(xx * radius, yy * radius, zz * radius); % 画一个球面做参考
set(h2,'edgecolor','none','facecolor','k','facealpha',1);
grid on
xlabel='x';
ylabel='y';
zlabel='z';
axis equal

%------------------------------------------------
function [points,b] = fitplane(points)
%x,y,z必须是列向量
% n=length(points);
x=points(:,1);
y=points(:,2);
z=points(:,3);
b = regress([z;-z], [[x;-x],[y;-y]]);
t=([x,y,z]*[b;-1])/(sum([b;-1].^2));
x=x-b(1)*t;
y=y-b(2)*t;
z=z+t;
points=[x,y,z];
end

%------------------------------------------------
% function [points] = fitplane(points)


