function NURBS = gen_Iee_Ine_surf(NURBS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate the Iee matrix, which relates element numbers and local edge number to the appropriate global egde numbers; %%%
%%% Generate the Ine matrix, which relates global egde number to the physical coordinates with dummy nodes. %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from NURBS
p = NURBS.p; uKnot = NURBS.uKnot; mcp = NURBS.mcp;
q = NURBS.q; vKnot = NURBS.vKnot; ncp = NURBS.ncp;
CP = NURBS.CP;

%% ===== Initialization =====
nel = (mcp - p) * (ncp - q); nedge_x = (mcp - p) * (ncp - q + 1); nedge_y = (ncp - q) * (mcp - p + 1); npoint = 101;
Iee = zeros(nel, 4); Ine = zeros(nedge_x + nedge_y, npoint, 3);

%% ===== Generate edge physical coordinates and edge indexes of elements ======
% --- Edge separations on v-direction (along u-direction) ---
for edge_y = 1:ncp-q+1  % Loop for number of edges on v-direction (vKnot ~= 0 && ~= 1)
    v = vKnot(edge_y+q); nj = find_knot_span(v,vKnot,ncp);
    for edge_x = 1:mcp-p  % Loop for number of elements on u-direction
        for i_point = 1:npoint
            u = uKnot(edge_x+p) + (uKnot(edge_x+p+1) - uKnot(edge_x+p))/(npoint - 1) * (i_point-1);
            ni = find_knot_span(u,uKnot,mcp);
            i_edge = (edge_y - 1)*(mcp-p) + edge_x;
            Ine(i_edge, i_point, :) = cal_phys_coord_surf(ni,p,u,uKnot,nj,q,v,vKnot,CP);
            
            if edge_y <= ncp-q
                i_e = (edge_y-1)*(ncp-q) + edge_x;
                Iee(i_e, 1) = i_edge; Iee(i_e, 3) = i_edge + (mcp-p);
            end
        end
    end
end

% --- Edge separations on u-direction (along v-direction) ---
for edge_x = 1:mcp-p+1  % Loop for number of edges along u-direction (uKnot ~= 0 && ~= 1)
    u = uKnot(edge_x+p); ni = find_knot_span(u,uKnot,mcp);    
    for edge_y = 1:ncp-q  % Loop for number of elements on v-direction
        for i_point = 1:npoint
            v = vKnot(edge_y+q) + (vKnot(edge_y+q+1) - vKnot(edge_y+q))/(npoint - 1) * (i_point-1);
            nj = find_knot_span(v,vKnot,ncp);
            i_edge = (edge_x - 1)*(ncp-q) + edge_y + nedge_x;
            Ine(i_edge, i_point, :) = cal_phys_coord_surf(ni,p,u,uKnot,nj,q,v,vKnot,CP);
            
            if edge_x <= mcp-p
                i_e = (edge_y-1)*(ncp-q) + edge_x;
                Iee(i_e, 2) = i_edge; Iee(i_e, 4) = i_edge + (ncp-q);
            end
        end
    end
end

%% ===== Save NURBS =====
NURBS.Iee = Iee; NURBS.Ine = Ine;
end
