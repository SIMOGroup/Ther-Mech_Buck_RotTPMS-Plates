function NURBS = gen_Ien_Inn_surf(NURBS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate the Ien matrix, which relates element numbers and local node numbers to the appropriate global node numbers; %%%
%%% Generate the Inn matrix, which relates global node number to the NURBS natural coordinates of the node. %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from NURBS
p = NURBS.p; mcp = NURBS.mcp; 
q = NURBS.q; ncp = NURBS.ncp;

%% ===== Initialization =====
nel = (mcp - p) * (ncp - q); nnode_per_el = (p +1) * (q +1); nnode = mcp * ncp;
Ien = zeros(nel, nnode_per_el); Inn = zeros(nnode, 2);

%% ===== Compute Ien and Inn =====
for j = 1:ncp  % Loop through control points in v direction
    for i = 1:mcp  % Loop through control points in u direction
        g = (j-1)*mcp + i;
        Inn(g,1) = i; Inn(g,2) = j;
        if i >= p+1 && j >= q+1
            i_e = (j-q-1)*(ncp-q) + i-p;
            for k_v = 1:q+1  % Loop through control points in v direction of one element
                for k_u = 1:p+1  % Loop through control points in v direction of one element
                    idx = (k_v-1)*(p+1) + k_u;
                    Ien(i_e,idx) = g - (k_v-1)*mcp - (k_u-1);
                end
            end
        end
    end
end

%% ===== Save NURBS =====
NURBS.Ien = Ien; NURBS.Inn = Inn; 
end
