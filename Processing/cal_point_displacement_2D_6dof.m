function [point_disp_norm] = cal_point_displacement_2D_6dof(IGA,Plate,phys_coord)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate deflection of specific point of rectangular plate with six-variable plate theory %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
uKnot = IGA.NURBS.uKnot; vKnot = IGA.NURBS.vKnot;
Ien = IGA.NURBS.Ien; Inn = IGA.NURBS.Inn;
Iee = IGA.NURBS.Iee; Ine = IGA.NURBS.Ine;
gcoord = IGA.NURBS.gcoord; nel = IGA.NURBS.nel;
ndof = IGA.params.ndof;
U = IGA.result.U;
norm_method = IGA.result.norm_method;

%% Used parameters from Plate
L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h;
shear_func = Plate.theory.shear_func; stretch_func = Plate.theory.stretch_func;
q_uniform = Plate.prob.q_uniform;

%% ===== Point deflection =====
phys_coord = [phys_coord(1), phys_coord(2), 0];
% --- Find element ---
for i_ele = 1:size(Iee,1)    
    % Search within element coordinate box
    edges = Iee(i_ele, :);
    box_min = reshape(min(min(Ine(edges,:,:),[],1),[],2),1,[]);
    box_max = reshape(max(max(Ine(edges,:,:),[],1),[],2),1,[]);
    
    for i_dim = 1:3
        if phys_coord(i_dim) >= box_min(i_dim) && phys_coord(i_dim) <= box_max(i_dim)
            check = true;
        else
            check = false;
            break
        end
    end
    if check
        c_ele = i_ele;
        break
    end
end

% --- Element parameters ---
sctr = Ien(c_ele,:);  % Control points indexes
ni = Inn(Ien(c_ele,1),1);  % Index of the element in parametric domain
nj = Inn(Ien(c_ele,1),2);
nn = length(sctr);  % Number of control points in the element
for idof = 1:ndof  % Dofs of control points
    sctrF(idof:ndof:ndof*(nn-1) + idof) = ndof.*(sctr-1) + idof;
end
nodes = gcoord(sctr,:);

% --- Point parameters ---
nat_coord = cal_nat_coord_surf(IGA.NURBS,ni,nj,nodes,phys_coord);  % Natural coordinate
nat_coord(1) = (uKnot(ni+1)-uKnot(ni))/2*nat_coord(1) + (uKnot(ni+1)+uKnot(ni))/2;  % Map the point to parametric domain
nat_coord(2) = (vKnot(nj+1)-vKnot(nj))/2*nat_coord(2) + (vKnot(nj+1)+vKnot(nj))/2;

% --- Kinematic matrices of NURBS shape function ---
[N, ~, ~, ~, ~, ~] = Kine_Shape_2nd_2D(IGA.NURBS,ni,nj,nat_coord(1),nat_coord(2));
                
% --- Stretching effect deformation function ---
[gz, ~] = compute_stretching_deformation_function(0,h,shear_func,stretch_func);

% --- Deflection ---
N0 = zeros(1,ndof*nn);
N0(1,3:ndof:ndof*nn) = N';
N0(1,6:ndof:ndof*nn) = gz*N';
cen_def = N0*U(sctrF);

% --- Normalization ---
switch norm_method
    case 1
        [E_m, nu_m, ~, ~, ~] = compute_temperature_dependent_material(Plate.mat.type, 300);
%         [E_m, nu_m, rho_m] = compute_material(Plate.mat.type);
        D_m = E_m*h^(3)/(12*(1-nu_m^2));
        norm = 100*D_m/(q_uniform*L^4);
        point_disp_norm = cen_def*norm;
    case 2
        norm = 1/h;
        point_disp_norm = cen_def*norm;
end
end
