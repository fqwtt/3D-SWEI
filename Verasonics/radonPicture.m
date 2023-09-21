function radonPicture(disImg,theta)
    figure();imagesc(disImg);hold on
    [R,xp]=radon(disImg,theta);
    [~,ind]=min(min(R));
    [~,i1]=min(R);
    i1=xp(i1);i1=i1(ind);
    loc= floor((size(disImg)+1)/2);
    x0=loc(2)+cosd(theta(ind))*i1;y0=loc(1)-sind(theta(ind))*i1;k=cotd(theta(ind));b=y0-k*x0;plot(1:90,k*(1:90)+b,'w');
    hold off
end

