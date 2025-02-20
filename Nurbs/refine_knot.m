function R = refine_knot(Knot, ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Refine 1D Knot vector by insert knot %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
R = [];
k = 1;
for i = 1:length(Knot)-1
    if Knot(i) ~= Knot(i+1)
        for j = 1:ref-1
            R(k) = j/ref * (Knot(i+1) - Knot(i)) + Knot(i);
            k = k + 1;
        end 
    end
end
end