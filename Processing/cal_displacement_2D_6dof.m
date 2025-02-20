function [NURBS_disp_norm] = cal_displacement_2D_6dof(IGA,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate deflection of all control points of rectangular plate with six-variable plate theory %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
CP = IGA.NURBS.CP;
ndof = IGA.params.ndof; 
U = IGA.result.U;
norm_method = IGA.result.norm_method;

%% Used parameters from Plate
L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h;
shear_func = Plate.theory.shear_func; stretch_func = Plate.theory.stretch_func;
q_uniform = Plate.prob.q_uniform;

%% ===== Initial deformation NURBS =====
NURBS_disp_norm = zeros(size(CP,1),size(CP,2));

%% ====== Plot deflection ======
for j = 1:size(CP,2)
    for i = 1:size(CP,1) 
        idx = (j-1)*size(CP,1) + i;
        
        % --- Stretching effect deformation function ---
        [gz, ~] = compute_stretching_deformation_function(0,h,shear_func,stretch_func);
        
        disp = U(ndof*(idx-1)+3) + gz*U(ndof*(idx-1)+6);
        
        % --- Normalization ---
        switch norm_method
            case 0 
                disp_norm = disp;
            case 1
                [E_m, nu_m, ~, ~, ~] = compute_temperature_dependent_material(Plate.mat.type, 300);
%                 [E_m, nu_m, ~, ~] = compute_basic_material(Plate.mat.type);
                D_m = E_m*h^(3)/(12*(1-nu_m^2));
                norm = 100*D_m/(q_uniform*L^4);
                disp_norm = disp*norm;
            case 2
                norm = 1/h;
                disp_norm = disp*norm;
        end
        NURBS_disp_norm(i,j) = disp_norm;
    end
end
end
