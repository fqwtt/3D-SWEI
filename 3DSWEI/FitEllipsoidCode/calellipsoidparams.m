function [Center,Axis,R] = calellipsoidparams(C)

M = [C(1) C(4)/2 C(5)/2;
     C(4)/2 C(2) C(6)/2;
     C(5)/2 C(6)/2 C(3)];
% 中心
Center = -0.5*[C(7) C(8) C(9)]*inv(M);

S = Center*M*Center'+1;
[U,V]=eig(M);

[~,indx] = max(abs(U(1,:)));
[~,indy] = max(abs(U(2,:)));
[~,indz] = max(abs(U(3,:)));

Axis = [sqrt(S/V(indx,indx)) sqrt(S/V(indy,indy)) sqrt(S/V(indz,indz))];
R = U(:,[indx,indy,indz]);
end
