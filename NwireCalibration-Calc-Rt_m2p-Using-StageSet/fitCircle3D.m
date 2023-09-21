function [C, n] = fitCircle3D(M)
[num, dim]=size(M);
 
L1=ones(num,1);
A=inv(M'*M)*M'*L1; % plane normal vector    
 
B=zeros((num-1)*num/2,3);
 
count=0;
for i=1:num-1
    for j=i+1:num   
        count=count+1;
        B(count,:)=M(j,:)-M(i,:);
    end    
end
 
L2=zeros((num-1)*num/2,1);
count=0;
for i=1:num-1
    for j=i+1:num
        count=count+1;
        L2(count)=(M(j,1)^2+M(j,2)^2+M(j,3)^2-M(i,1)^2-M(i,2)^2-M(i,3)^2)/2;
    end
end
 
D=zeros(4,4);
D(1:3,1:3)=(B'*B);
D(4,1:3)=A';
D(1:3,4)=A;
 
L3=[B'*L2;1];

%C = pinv(D')*(L3);
%C = pinv(D') * L3 + (eye(4) - pinv(D') * D');
%C = lsqminnorm(D',L3);
%C=inv(D')*(L3);  % circle center 
C = D' \ L3;

C=C(1:3);
 
radius=0;
for i=1:num
    tmp=M(i,:)-C';
    radius=radius+sqrt(tmp(1)^2+tmp(2)^2+tmp(3)^2);
end
r=radius/num;            %  radius
 
n = A / norm(A);

%% Drawing
%figure;

% Draw data points
h1=plot3(M(:,1),M(:,2),M(:,3),'*'); 

% Draw whole circle
hold on;
c = C;
theta=(0:2*pi/100:2*pi)';   
a=cross(n,[1 0 0]);          
if ~any(a)                  
    a=cross(n,[0 1 0]);
end
b=cross(n,a);     
a=a/norm(a);      
b=b/norm(b);      
 
c1=c(1)*ones(size(theta,1),1);
c2=c(2)*ones(size(theta,1),1);
c3=c(3)*ones(size(theta,1),1);
 
x=c1+r*a(1)*cos(theta)+r*b(1)*sin(theta);  
y=c2+r*a(2)*cos(theta)+r*b(2)*sin(theta);  
z=c3+r*a(3)*cos(theta)+r*b(3)*sin(theta);  
 
h2=plot3(x,y,z,'-r');

% Draw fitted circle center
plot3(c(1), c(2), c(3), 'b*', 'MarkerSize', 1);

axis equal;
