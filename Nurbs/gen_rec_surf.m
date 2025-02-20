function NURBS = gen_rec_surf(Plate, NURBS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate NURBS surf for 2D rectangle %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
deg = NURBS.deg; ref = NURBS.ref; 

%% Used parameters from Plate
L = Plate.geo.L; W = Plate.geo.W;

%% ===== Rectangle mesh generation ======
% --- Coarse mesh ---
p = 1; q = 1; uKnot = [0, 0, 1, 1]; vKnot = [0, 0, 1, 1];
CP(:,:,1) = [0, 0; L, L]; CP(:,:,2) = [0, W; 0, W]; 
CP(:,:,3) = [0, 0; 0, 0]; CP(:,:,4) = [1, 1; 1, 1];

% --- k-Refinement ---
[CP, uKnot, vKnot, p, q] = refine_p_surf(p,q,uKnot,vKnot,CP,deg-p,deg-q);  % p-Refinement
[CP, uKnot, vKnot] = refine_h_surf(p,q,uKnot,vKnot,CP,ref,ref);  % h-Refinement

%% ====== Reduced to 2D mesh (xy plane) ======
mcp = length(CP(:,1,1)); ncp = length(CP(1,:,1));
B_net(:,:,1)= CP(:,:,1); B_net(:,:,2)= CP(:,:,2); B_net(:,:,3)= CP(:,:,4);

%% ===== Global Coordinate of Control Points ======
gcoord(:,1) = reshape(B_net(:,:,1),mcp*ncp,1);
gcoord(:,2) = reshape(B_net(:,:,2),mcp*ncp,1);
gcoord(:,3) = reshape(B_net(:,:,3),mcp*ncp,1);

%% ===== Save NURBS =====
NURBS.p = p; NURBS.uKnot = uKnot; NURBS.mcp = mcp;
NURBS.q = q; NURBS.vKnot = vKnot; NURBS.ncp = ncp;
NURBS.CP = CP; NURBS.B_net = B_net; NURBS.gcoord = gcoord;
NURBS.deg = deg; NURBS.ref = ref; 

end
