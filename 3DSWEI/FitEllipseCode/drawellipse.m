function [x,y] = drawellipse(W)

% 获取椭圆参数
[Center,Axis,Theta] = calellipseparams(W);

% 绘制椭圆
funs = @(x,y) W(1)*x.^2 + W(2)*x.*y + W(3)*y.^2 + W(4)*x + W(5)*y + W(6);
h=fimplicit(funs,'LineWidth',2);
x=get(h,'Xdata')';
y=get(h,'Ydata')';
% 绘制长短轴
Majcoor = [-Axis(1),0; Axis(1),0];
Mincoor = [0,-Axis(2); 0,Axis(2)];
RM = [cos(Theta),-sin(Theta);sin(Theta),cos(Theta)];
Majcoor = Majcoor*RM'+Center;
Mincoor = Mincoor*RM'+Center;

line(Majcoor(:,1),Majcoor(:,2),'Color','r','LineWidth',3)
line(Mincoor(:,1),Mincoor(:,2),'Color','g','LineWidth',3)
plot(Center(1),Center(2),'y.','MarkerSize',15)

end




