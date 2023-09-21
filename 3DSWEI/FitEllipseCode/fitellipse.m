function W = fitellipse(x,y)

% 构造矩阵
D = [x.*x, x.*y, y.*y, x, y,ones(size(x))];
S = D'*D;
G = zeros(6);
G(1,3) = 2; G(3,1) = 2; G(2,2) = -1;

% 求解
[vec, val] = eig(S\G);
[~, idx] = find(val>0&~isinf(val));
W = vec(:,idx);
W = sqrt(1/(W'*S*W))*W;
end

