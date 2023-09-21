%% 导入数据
data_w2c=importdata('..\1123\2\fileRtBack.txt');
data_t2w=importdata('..\..\PNN\PNN\fileRt_t2w_L7_3.txt');
n=length(data_w2c)/2;
n=630;

%% 调整数据格式
R_w2c=zeros(3,3,n);
t_w2c=zeros(3,1,n);
for i=1:n
    R_w2c(:,:,i)=reshape(data_w2c(2*i-1,1:9),3,3)';
    t_w2c(:,:,i)=data_w2c(2*i,1:3)';
end
R_t2w=reshape(data_t2w(1,1:9),3,3)';
t_t2w=data_t2w(2,1:3)';



tran1=[[0,0,-1];[-1,0,0];[0,1,0]];
R_c2c1=zeros(3,3,n);

tran2=[[1,0,0];[0,0,0];[0,1,-1]];
R_c2c2=zeros(3,3,n);

for i=1:n
    R_c2c1(:,:,i)=tran1*R_w2c(:,:,i);
    R_c2c2(:,:,i)=tran2*R_w2c(:,:,i);
end
%% 计算
x1=zeros(1,n);
y1=zeros(1,n);
z1=zeros(1,n);
x2=zeros(1,n);
y2=zeros(1,n);
z2=zeros(1,n);
x3=zeros(1,n);
y3=zeros(1,n);
z3=zeros(1,n);

newks=min(abs(ks(:,1)),abs(ks(:,2)));
for i=1:n
    tmp=R_w2c(:,:,i)*R_t2w;
    x1(i)=tmp(1,1)*newks(i);
    y1(i)=tmp(2,1)*newks(i);
    z1(i)=tmp(3,1)*newks(i);
  
%     tmp=R_c2c1(:,:,i)*R_t2w;
%     x2(i)=-tmp(1,1)*new_ks(i,1);
%     y2(i)=-tmp(2,1)*new_ks(i,1);
%     z2(i)=-tmp(3,1)*new_ks(i,1);
% 
%     tmp=R_c2c2(:,:,i)*R_t2w;
%     x3(i)=tmp(1,1)*new_ks(i,1);
%     y3(i)=tmp(2,1)*new_ks(i,1);
%     z3(i)=tmp(3,1)*new_ks(i,1);
end

%% 绘图
X=-3:0.1:3;
[X,Y]=meshgrid(X);
plot3(x1,y1,z1,'*');
xlabel='x';
ylabel='y';
zlabel='z';
axis equal
grid;