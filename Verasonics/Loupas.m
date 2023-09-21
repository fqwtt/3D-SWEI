function [DisData] = Loupas(IBuffer, QBuffer)

    %% 数据预处理
    c=evalin('base','c');
    fc=evalin('base','fc');
    na = evalin('base','na');
    PData = evalin('base','PData');
    Depth = evalin('base','Depth');
    Width = evalin('base','Width');
    ROILA = PData(2).Region.PixelsLA(1);
    [h,~,~]=size(IBuffer);
    point=[mod(ROILA,h),floor(ROILA/h)];
    IBuffer=IBuffer(point(1)+1:point(1)+Depth,point(2)+1:point(2)+Width,:);
    QBuffer=QBuffer(point(1)+1:point(1)+Depth,point(2)+1:point(2)+Width,:);

    %% 可变参数
    M=33;
    N=7;
    kernel2D = ...
        [   0.0073    0.0208    0.0294    0.0208    0.0073;
            0.0208    0.0589    0.0833    0.0589    0.0208;
            0.0294    0.0833    0.1179    0.0833    0.0294;
            0.0208    0.0589    0.0833    0.0589    0.0208;
            0.0073    0.0208    0.0294    0.0208    0.0073;];
    %% 计算部分
    Im1=QBuffer(:,:,1:na-1).*IBuffer(:,:,2:na) - IBuffer(:,:,1:na-1).*QBuffer(:,:,2:na);
    Re1=IBuffer(:,:,1:na-1).*IBuffer(:,:,2:na) + QBuffer(:,:,1:na-1).*QBuffer(:,:,2:na);
    Im2=QBuffer(1:Depth-1,:,:).*IBuffer(2:Depth,:,:) - IBuffer(1:Depth-1,:,:).*QBuffer(2:Depth,:,:);
    Re2=IBuffer(1:Depth-1,:,:).*IBuffer(2:Depth,:,:) + QBuffer(1:Depth-1,:,:).*QBuffer(2:Depth,:,:);
    kerm1=ones(M,1);kerm2=ones(M-1,1);
    kern1=ones(1,1,N);kern2=ones(1,1,N-1);
    Im1=convn(convn(Im1,kern2,'valid'),kerm1,'valid');
    Re1=convn(convn(Re1,kern2,'valid'),kerm1,'valid');
    Im2=convn(convn(Im2,kern1,'valid'),kerm2,'valid');
    Re2=convn(convn(Re2,kern1,'valid'),kerm2,'valid');
    
    DisData=convn(atan2(Im1 , Re1) ./ (1 + atan2(Im2 , Re2)/(2 * pi)) *c / (4 * pi * fc),kernel2D,'same');
%     DisData=atan2(Im1 , Re1) ./ (1 + atan2(Im2 , Re2)/(2 * pi)) *c / (4 * pi * fc);
end
    
    
    
    