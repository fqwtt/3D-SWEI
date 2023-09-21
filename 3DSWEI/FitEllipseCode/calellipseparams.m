function [Center,Axis,Theta] = calellipseparams(W)

a = W(1); 
b = W(2); 
c = W(3); 
d = W(4); 
e = W(5); 
f = W(6);

% 中心
cx = (b*e-2*c*d)/(4*a*c-b^2);
cy = (b*d-2*a*e)/(4*a*c-b^2);
Center = [cx,cy];

% 长短轴  
MA1 = sqrt(2*(a*cx^2+c*cy^2+b*cx*cy-f)/(a+c+sqrt((a-c)^2+b^2)));
MA2= sqrt(2*(a*cx^2+c*cy^2+b*cx*cy-f)/(a+c-sqrt((a-c)^2+b^2)));

Axis = [max(MA1,MA2),min(MA1,MA2)];

% 长轴倾角
if b==0
    if f*a>f*c
        Theta = 0;
    else  
        Theta = pi/2;
    end
else
    if f*a>f*c
        alpha = atan((a-c)/b);
        if alpha<0
            Theta = 0.5*(-pi/2-alpha);
        else
            Theta = 0.5*(pi/2-alpha);
        end
    else
        alpha = atan((a-c)/b);
        if alpha<0
            Theta = pi/2+0.5*(-pi/2-alpha);
        else
            Theta = pi/2+0.5*(pi/2-alpha);
        end
    end
end

end
