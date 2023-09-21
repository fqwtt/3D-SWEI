function Yang = radonSum(dis,push)
    Ledge=evalin('base','Ledge');
    Redge=evalin('base','Redge');
    unitTime=evalin('base','unitTime');
    unitDistance=evalin('base','unitDistance');
    nf=evalin('base','nframe');
    [depth,width,~]=size(dis);
    Yang=zeros(depth,width);
    aver_ker=ones(8,1)/8; %在轴向上求平均
    dis=convn(dis,aver_ker,'same');
    thetaR=[10:80];thetaL=[100:170];
    for n=1:length(push)
        Yang=zeros(depth,width);
        for i=1:depth
            for j=[Ledge:push(n).focus]%push(n).LRange
                disImg=squeeze(dis(i,j-15:j+5,1+nf*(n-1):90+nf*(n-1)));
                [R,~]=radon(disImg,thetaL);
                [~,ind]=min(min(R));
                Yang(i,j)=thetaL(ind);
%                 disImg=squeeze(dis(200,36:164,1+nf*(n-1):90+nf*(n-1)));radonPicture(disImg,thetaL);
            end
            for j=[push(n).focus+1:Redge]%push(n).RRange
                disImg=squeeze(dis(i,j-5:j+15,1+nf*(n-1):90+nf*(n-1)));
                [R,~]=radon(disImg,thetaR);
                [~,ind]=min(min(R));
                Yang(i,j)=thetaR(ind);
%                 disImg=squeeze(dis(200,36:164,1+nf*(n-1):90+nf*(n-1)));radonPicture(disImg,thetaR);
            end
        end
        Yang=(cotd(Yang)*unitDistance/unitTime).^2*3000;
        figure()
        imagesc(Yang(:,37:164),[0,10000]);
        title('Image of elasticity modulus');
        xlabel('lateral position(wavelength)');ylabel('axial position(wavelength)');
        set(gca,'XTick',linspace(1,128,7),'XTickLabel',round(linspace(-32,32,7)),...
                'YTick',linspace(1,400,11),'YTickLabel',round(linspace(0,100,11)));
    end
%     Yang=(cotd(Yang)*unitDistance/unitTime).^2*3000;
end

