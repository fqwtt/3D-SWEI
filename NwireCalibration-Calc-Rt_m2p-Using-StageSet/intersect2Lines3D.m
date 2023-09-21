function K_Pnt = intersect2Lines3D(v1, p1, v2, p2)
L1_Dir = v1;

L1_Dir = L1_Dir./norm(L1_Dir);

P1_Pnt = p1;

L2_Dir = v2;

L2_Dir = L2_Dir./norm(L2_Dir);

P2_Pnt = p2;

P1P2_Vec = P2_Pnt-P1_Pnt;

P2P3_Norm = norm(cross(P1P2_Vec,L1_Dir));

P3_Pnt = P1_Pnt + (P1P2_Vec*L1_Dir').*L1_Dir;

CosTheta= abs(L1_Dir*L2_Dir');
K_Pnt = [];
if CosTheta < 1e-7
    
    K_Pnt = P3_Pnt;
end

if CosTheta > 1e-7
    TanTheta = (1-CosTheta^2)^0.5/CosTheta;
    KP3_Norm = P2P3_Norm/TanTheta;
    
    K1_Pnt = P3_Pnt + KP3_Norm.*L1_Dir;
    K2_Pnt = P3_Pnt - KP3_Norm.*L1_Dir;
    
    P2K1_Vec = K1_Pnt - P2_Pnt;
    P2K2_Vec = K2_Pnt - P2_Pnt;
    
    D1 = norm(cross(P2K1_Vec,L2_Dir));
    D2 = norm(cross(P2K2_Vec,L2_Dir));
    
    if D1 < D2
        K_Pnt = K1_Pnt;
    else
        K_Pnt = K2_Pnt;
    end
    
end

end