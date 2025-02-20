function [phys_coord] = cal_phys_coord_surf(ni,p,u,uKnot,nj,q,v,vKnot,CP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate physical coordinates of NURBS surf of a natural knot point %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
Nu = eval_basis_func(ni,p,u,uKnot); 
Nv = eval_basis_func(nj,q,v,vKnot);

SumNw = 0;
for j = 0:q
    for i = 0:p
        SumNw = Nu(i+1)*Nv(j+1)*CP(ni-p+i,nj-q+j,4) + SumNw;
    end
end

phys_coord(1:3) = 0;
for j = 0:q 
    for i = 0:p
        phys_coord(1) = Nu(i+1)*Nv(j+1)*CP(ni-p+i,nj-q+j,4)*CP(ni-p+i,nj-q+j,1)/SumNw + phys_coord(1);
        phys_coord(2) = Nu(i+1)*Nv(j+1)*CP(ni-p+i,nj-q+j,4)*CP(ni-p+i,nj-q+j,2)/SumNw + phys_coord(2);
        phys_coord(3) = Nu(i+1)*Nv(j+1)*CP(ni-p+i,nj-q+j,4)*CP(ni-p+i,nj-q+j,3)/SumNw + phys_coord(3); 
    end
end
end
