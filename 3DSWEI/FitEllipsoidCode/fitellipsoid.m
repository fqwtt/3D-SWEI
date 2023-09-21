function C = fitellipsoid(x,y,z)

% 构造矩阵
X = [x.*x y.*y z.*z x.*y x.*z y.*z x y z];
Y = ones(length(x),1);

% 求解
C = inv(X'*X)*X'*Y;
end
