function p = polyfit_wt(x,y)
    x = x(:);
    y = y(:);
    V(:,2) = ones(length(x),1);
    V(:,1) = x.*V(:,2);
    [Q,R] = qr(V);
    p = R\(Q'*y);
end

