function NURBS = gen_FE_approx_surf(NURBS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Finite element approximated surf %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
p = NURBS.p; uKnot = NURBS.uKnot; mcp = NURBS.mcp;
q = NURBS.q; vKnot = NURBS.vKnot; ncp = NURBS.ncp;
CP = NURBS.CP;

%% ===== Initialization =====
npoint = 101;
FE_grid = zeros(npoint, npoint, 3);

%% ===== Generate edge physical coordinates and edge indexes of elements ======
for i_point_y = 1:npoint
    v = vKnot(q+1) + (vKnot(ncp+1) - vKnot(q+1))/(npoint - 1) * (i_point_y-1);
    nj = find_knot_span(v,vKnot,ncp);
    for i_point_x = 1:npoint
        u = uKnot(p+1) + (uKnot(mcp+1) - uKnot(p+1))/(npoint - 1) * (i_point_x-1);
        ni = find_knot_span(u,uKnot,mcp);
        FE_grid(i_point_x, i_point_y, :) = cal_phys_coord_surf(ni,p,u,uKnot,nj,q,v,vKnot,CP);
    end
end

%% ===== Save NURBS =====
NURBS.FE_grid = FE_grid;
