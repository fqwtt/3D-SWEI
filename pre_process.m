%% 参数
tic
nline        = 1;
nframe       = 100;
Depth        = 400;
Width        = 200;
thread       = 50;
c            = 1480;
fc           = 5.208;
Ledge        = 36;
Redge        = 164;
unitTime     = 0.1;
unitDistance = c / (fc * 1000) / 2;  %0.1478
na           = nline * nframe;
model        =[0,0];
%% 获得IQ数据
idata=IData(2);qdata=QData(2);
idata=squeeze(idata{1});qdata=squeeze(qdata{1});

%% 使用Kasai/Loupas算法获得位移数据
% dis=Kasai2(idata,qdata);
dis=Loupas(idata,qdata);
% dis=Loupas(IData,QData);


%% 低通滤波
% dis=-dis./min(dis(:,:,:),[],[2]);
% dis(dis<-1)=0;
dis=lowpass2(dis,thread);

%% 轴向平均处理
dis=convn(dis,ones(50,1)/50,'same');

%% 使用三次样条插值将时间分辨率提高5倍
% dis = spline([1:na-2], dis, [1:0.2:na-2]);

%% 使用三次样条插值将横向分辨率提高5倍
% dis=permute(spline([1:Width],...
%                    permute(dis,[1,3,2]),...
%                    [1:0.2:Width]),[1,3,2]);

%% Radon和变换求速度

push(1).focus=40; push(1).LRange=[36:50];push(1).RRange=[70:100];
push(2).focus=160;push(2).LRange=[50:90];push(2).RRange=[110:150];
% push(3).focus=140;push(3).LRange=[100:130];push(3).RRange=[150:164];
% push(4).focus=150;push(4).left=[108:120];push(4).right=[160:164];
% dis(:,1:36,:)=0;
% dis(:,164:200,:)=0;
% dis=-dis./min(dis(:,36:164,:),[],[1,2]); %归一化
% tic
% Yang=radonSum(dis,push);
% toc

%% Ransac求速度

Yang=Ransac_wt(dis);
toc



