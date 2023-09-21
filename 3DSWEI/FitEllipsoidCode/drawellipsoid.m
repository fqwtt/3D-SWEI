function drawellipsoid(C)

% 获取椭圆参数
[Center,Axis,R] = calellipsoidparams(C);
xc = Center(1);
yc = Center(2);
zc = Center(3);
xr = Axis(1);
yr = Axis(2);
zr = Axis(3);
% 绘制椭圆
pointNum = 500;
[x,y,z] = ellipsoid(xc,yc,zc,xr,yr,zr,pointNum);
temp = R*[x(:),y(:),z(:)]';
xx = reshape(temp(1,:),[pointNum+1,pointNum+1]);
yy = reshape(temp(2,:),[pointNum+1,pointNum+1]);
zz = reshape(temp(3,:),[pointNum+1,pointNum+1]);
cc = sqrt(xx.*xx + yy.*yy + zz.*zz);
surf(xx,yy,zz,cc,'EdgeColor','none','FaceAlpha',1.0);
axis equal;
end




