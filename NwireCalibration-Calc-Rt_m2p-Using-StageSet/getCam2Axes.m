function [R_c2a, t_c2a] = getCam2Axes(py, diry, pz, dirz, shifty, isYaccurate)

diry = reshape(diry, [1,3]);
dirz = reshape(dirz, [1,3]);

% ACCRATEY is true -- axis Y is considered accurate, and axis Z is
% translated to meet axis Y; if set to false, axis Z will be considered
% accurate.
ACCURATEY = isYaccurate;

if ACCURATEY == true
    vecy = diry / norm(diry);
    vecz0 = dirz / norm(dirz);
    vecx = cross(vecy, vecz0);
    vecx = vecx / norm(vecx);
    vecz = cross(vecx, vecy);
    vecz = vecz / norm(vecz);
    
    syms m;
    p = pz + m * vecx; % translate axis z to meet axis y (p is not their meeting point!)
    L = [(py - p); vecz0; vecy]; % means the 3 vectors are coplanar
    func = det(L);
    m = solve(func);
    m = double(m);
    p = double(subs(p));
    % the original point is the intersecting point of the two axes (after
    % translation)
    p_orig = intersect2Lines3D(vecz0, p, vecy, py);
    
else
    vecy0 = diry / norm(diry);
    vecz = dirz / norm(dirz);
    vecx = cross(vecy0, vecz);
    vecx = vecx / norm(vecx);
    vecy = cross(vecz, vecx);
    vecy = vecy / norm(vecy);
    
    syms m;
    p = py + m * vecx; % translate axis y to meet axis z (p is not their meeting point!)
    L = [(py - p); vecy0; vecz]; % means the 3 vectors are coplanar
    func = det(L);
    m = solve(func);
    m = double(m);
    p = double(subs(p));
    % the original point is the intersecting point of the two axes (after
    % translation)
    p_orig = intersect2Lines3D(vecz, pz, vecy0, p);
end

%p_orig2  = p_orig + vecy * shifty;
p_orig2  = p_orig + vecz * shifty;

R_c2a = [vecx; vecy; vecz];
t_c2a = -p_orig2';

%% Figure the translation result

%figure;
t = -100 : 0.5 : 100;

if ACCURATEY == true 
    p_z0 = p_orig - m * vecx;
    for i = 1: 3
        axY(:,i) = p_orig(i) + t .* vecy(i);
        axZ1(:,i) = p_orig(i) + t .* vecz0(i);
        axZ0(:,i) = p_z0(i) + t .* vecz0(i);
    end
    
    plot3(axY(:,1),axY(:,2), axY(:,3),'k');
    hold on;
    plot3(axZ0(:,1),axZ0(:,2), axZ0(:,3),'r');
    plot3(axZ1(:,1),axZ1(:,2), axZ1(:,3),'--r');

else
    p_y0 = p_orig - m * vecx;
    for i = 1: 3
        axY0(:,i) = p_y0(i) + t .* vecy0(i);
        axY1(:,i) = p_orig(i) + t .* vecy0(i);
        axZ(:,i) = p_orig(i) + t .* vecz(i);
    end
    
    plot3(axZ(:,1),axZ(:,2), axZ(:,3),'k');
    hold on;
    plot3(axY0(:,1),axY0(:,2), axY0(:,3),'r');
    plot3(axY1(:,1),axY1(:,2), axY1(:,3),'--r');
end

plot3(p_orig(1), p_orig(2), p_orig(3),'.')
axis equal;

end