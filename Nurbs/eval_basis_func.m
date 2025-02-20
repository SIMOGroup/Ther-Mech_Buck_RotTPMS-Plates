function N = eval_basis_func(i,p,u,Knot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Evaluate BSpline basis finction at a natural knot point in 1D Knot vector %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Source: Algorithm from Piegl, Les. "The NURBS Book". Springer-Verlag: Berlin 1995; pp. 72-73.
left = zeros(1,p+1); right = zeros(1,p+1);
ndu = zeros(p+1,p+1); N = zeros(1,p+1);

ndu(1,1) = 1;
for j = 1:p
    left(j+1) = u - Knot(i+1-j);
    right(j+1) = Knot(i+j) - u;
    saved = 0;
    for r = 0:j-1
        ndu(j+1,r+1) = right(r+2) + left(j-r+1);
        temp = ndu(r+1,j) ./ ndu(j+1,r+1);
        ndu(r+1,j+1) = saved + right(r+2) .* temp;
        saved = left(j-r+1) .* temp;
    end 
    ndu(j+1,j+1) = saved;
end 

for j = 0:p
    N(1,j+1) = ndu(j+1,p+1);
end

return
