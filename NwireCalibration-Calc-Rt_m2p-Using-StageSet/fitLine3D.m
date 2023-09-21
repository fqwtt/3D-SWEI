function [xyz0, direction] = fitLine3D(x, y, z, isFigure)
xyz0(1)=mean(x);
xyz0(2)=mean(y);
xyz0(3)=mean(z);

lineData = [x, y, z];
centeredLine=bsxfun(@minus,lineData,xyz0);
[U,S,V]=svd(centeredLine);
direction=V(:,1);

%% Drawing
if isFigure
    figure;
    
    plot3(x, y, z,'*');
    hold on;
    
    t=-100:0.1:100;
    xx=xyz0(1)+direction(1)*t;
    yy=xyz0(2)+direction(2)*t;
    zz=xyz0(3)+direction(3)*t;
    plot3(xx,yy,zz);
    axis equal;
end
