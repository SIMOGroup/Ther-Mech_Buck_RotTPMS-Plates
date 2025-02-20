function [Kg] = cal_Geometric_Stiffness_Matrices_2D_5dof(IGA,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate geometric stiffness matrices with five-variable plate theory %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
uKnot = IGA.NURBS.uKnot; vKnot = IGA.NURBS.vKnot;
Inn = IGA.NURBS.Inn; Ien = IGA.NURBS.Ien; nel = IGA.NURBS.nel;
sdof = IGA.params.sdof; ndof = IGA.params.ndof; nGauss = IGA.params.nGauss; 

%% Used parameters from Plate
buck_case = Plate.prob.buck_case;

%% ===== Initial stiffness matrices =====
Kg = sparse(sdof,sdof);
switch buck_case
    case 1  % X
        N0 = [1 0;0 0];  % Initial stress matrix
    case 2  % Y
        N0 = [0 0;0 1];
    case 3  % X+Y
        N0 = [1 0;0 1];
end

%% ===== Initial Gauss integration =====
[Pg, Wg] = gen_Gauss_point(1,-1,nGauss);
nel_nza = 0; % Number of non-zero-area elements
tol = 1e-8;

%% ===== Stiffness matrices =====
for iel = 1:nel  % Loop over the elements
    % --- Element parameters ---
    sctr = Ien(iel,:);  % Control points indexes
    nn = length(sctr);  % Number of control points in the element
    for idof = 1:ndof  % Dofs of control points
        sctrK(idof:ndof:ndof*(nn-1) + idof) = ndof.*(sctr-1) + idof;
    end
    ni = Inn(Ien(iel,1),1);  % Index of the element in parametric domain
    nj = Inn(Ien(iel,1),2);
    
    % --- Gauss integration ---
    if abs(uKnot(ni)-uKnot(ni+1)) > tol && abs(vKnot(nj)-vKnot(nj+1)) > tol  % Check if the current element has nonzero area in the parametric domain
        nel_nza = nel_nza + 1;
        detJ2_xi = (uKnot(ni+1) - uKnot(ni))/2;  % Mapping parametric domain into natural domain of [[-1, 1]; [-1, 1]]
        detJ2_eta = (vKnot(nj+1) - vKnot(nj))/2;

        for iGauss = 1: nGauss  % Loop over the integration points
            for jGauss = 1: nGauss
                % Gauss points & weights
                gpt_xi = Pg(iGauss); gwt_xi = Wg(iGauss);
                gpt_eta = Pg(jGauss); gwt_eta = Wg(jGauss);
                
                % Map the point to parametric domain
                gpt_xi = (uKnot(ni+1)-uKnot(ni))/2*gpt_xi + (uKnot(ni+1)+uKnot(ni))/2; gwt_xi = gwt_xi*detJ2_xi;
                gpt_eta = (vKnot(nj+1)-vKnot(nj))/2*gpt_eta + (vKnot(nj+1)+vKnot(nj))/2; gwt_eta = gwt_eta*detJ2_eta;
                gwt = gwt_xi*gwt_eta;
                
                % Kinematic matrices of NURBS shape function
                [N, dNdxi, dNdxi2, dNdxy, dN2dxy, detJ1] = Kine_Shape_2nd_2D(IGA.NURBS,ni,nj,gpt_xi,gpt_eta);
                
                % Stiffness matrix
                Bg = zeros(2,ndof*nn);  % Geometric
                Bg(1,3:ndof:ndof*nn) = dNdxy(:,1)';
                Bg(2,3:ndof:ndof*nn) = dNdxy(:,2)';
                
                kg = Bg'*N0*Bg*gwt*detJ1;
                
                Kg(sctrK,sctrK) = Kg(sctrK,sctrK) + kg;
            end
        end
    end
end
end
