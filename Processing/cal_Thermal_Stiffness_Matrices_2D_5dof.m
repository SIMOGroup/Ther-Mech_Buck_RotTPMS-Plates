function [K_th] = cal_Thermal_Stiffness_Matrices_2D_5dof(IGA,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate thermal stiffness matrices with five-variable plate theory %%%
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
S_th = Plate.mat_mat.S_th;

%% ===== Initial stiffness matrices =====
K_th = sparse(sdof,sdof);

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
                
                % Geometric stiffness matrix
                Bgx0 = zeros(3,ndof*nn);  % (u_0, v_0, w_0)
                Bgx0(1,1:ndof:ndof*nn) = dNdxy(:,1)';
                Bgx0(2,2:ndof:ndof*nn) = dNdxy(:,1)';
                Bgx0(3,3:ndof:ndof*nn) = dNdxy(:,1)';
                
                Bgx1 = zeros(3,ndof*nn);  % (-w_{0,x}, -w_{0,y}, 0)
                Bgx1(1,3:ndof:ndof*nn) = -dN2dxy(:,1)';
                Bgx1(2,3:ndof:ndof*nn) = -dN2dxy(:,3)';
                
                Bgx2 = zeros(3,ndof*nn);  % (\beta_{x}, \beta_{y}, 0)
                Bgx2(1,4:ndof:ndof*nn) = dNdxy(:,1)';
                Bgx2(2,5:ndof:ndof*nn) = dNdxy(:,1)';
                
                Bgy0 = zeros(3,ndof*nn);  % (u_0, v_0, w_0)
                Bgy0(1,1:ndof:ndof*nn) = dNdxy(:,2)';
                Bgy0(2,2:ndof:ndof*nn) = dNdxy(:,2)';
                Bgy0(3,3:ndof:ndof*nn) = dNdxy(:,2)';
                
                Bgy1 = zeros(3,ndof*nn);  % (-w_{0,x}, -w_{0,y}, 0)
                Bgy1(1,3:ndof:ndof*nn) = -dN2dxy(:,3)';
                Bgy1(2,3:ndof:ndof*nn) = -dN2dxy(:,2)';
                
                Bgy2 = zeros(3,ndof*nn);  % (\beta_{x}, \beta_{y}, 0)
                Bgy2(1,4:ndof:ndof*nn) = dNdxy(:,2)';
                Bgy2(2,5:ndof:ndof*nn) = dNdxy(:,2)';
                
                Bgx{1} = Bgx0; Bgx{2} = Bgx1; Bgx{3} = Bgx2;
                Bgy{1} = Bgy0; Bgy{2} = Bgy1; Bgy{3} = Bgy2;
                k_th = 0;
                for i = 1:3
                    for j = 1:3
                        k_th = k_th + S_th{i,j}(1,1)*Bgx{i}'*Bgx{j}*gwt*detJ1;
                        k_th = k_th + S_th{i,j}(2,2)*Bgy{i}'*Bgy{j}*gwt*detJ1;
                        k_th = k_th + 2 * S_th{i,j}(1,2)*Bgx{i}'*Bgy{j}*gwt*detJ1;
                    end
                end
                
                K_th(sctrK,sctrK) = K_th(sctrK,sctrK) + k_th;
            end
        end
    end
end
end
