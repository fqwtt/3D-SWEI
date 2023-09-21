sampleSize = 2; % number of points to sample per trial
maxDistance = 5; % max allowable distance for inliers
fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
evalLineFcn = @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);   % distance evaluation function

f=figure;
idx=8;
L=106;
R=160;
L2=41;
R2=95;
filepath = 'E:\nju307_wt\SWI_m\1123\2\test\';
filenumber = length(dir(strcat(filepath,'*.mat')));
% filenumber = 50;

% flag用来选择加载/显示数据
flag=1;

if flag
    data=zeros(200,44,filenumber); %#ok<UNRCH> 
    for i=1:filenumber
        load(strcat(filepath,string(i),'.mat'));
        data(:,:,i) = res1;
    end
end
ks=zeros(filenumber,4);
points=zeros(R-L+1,2);
for i =1:filenumber
    disp(i);
    subplot(1,3,1)
    temp=data(:,:,i);
    [~,b]=max(temp(:,:),[],2);
    plot(1:length(b),b);
    hold on;
    % 画左图右边的红线
    [~,b]=max(temp(L:R,:),[],2);
    points(:,1)=L:R;
    points(:,2)=b;
    [modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn,sampleSize,maxDistance);
    ks(i,2)=1/modelRANSAC(1);
    inlierPts = points(inlierIdx,:);
    x = [min(inlierPts(:,1)) max(inlierPts(:,1))];
    y = modelRANSAC(1)*x + modelRANSAC(2);
    plot(x, y, 'r-')
    hold on;

    % 画左图左边的红线
    [~,b]=max(temp(L2:R2,:),[],2);
    points(:,1)=L2:R2;
    points(:,2)=b;
    [modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn,sampleSize,maxDistance);
    ks(i,1)=1/modelRANSAC(1);
    inlierPts = points(inlierIdx,:);
    x1 = [min(inlierPts(:,1)) max(inlierPts(:,1))];
    y1 = modelRANSAC(1)*x1 + modelRANSAC(2);
    plot(x1, y1, 'r-')
    hold off;
    
    %imagesc
    subplot(1,3,2)
%     imagesc(temp);hold on;
    radonPicture(temp);hold on;
    plot(y,x,'w-');hold on;
    plot(y1,x1,'w-');hold off;
    


    % 画radon图
    subplot(1,3,3)
    modelRandon=radonPicture(temp(50:164,1:40));hold on;
    ks(i,3)=modelRandon;

%     waitforbuttonpress;
    if f.CurrentCharacter=='a'
        ks(i,4)=ks(i,1);
    end

    if f.CurrentCharacter=='d'
        ks(i,4)=ks(i,2);
    end

    if f.CurrentCharacter=='s'
        ks(i,4)=ks(i,3);
    end

    if f.CurrentCharacter=='q'
        break
    end
end


