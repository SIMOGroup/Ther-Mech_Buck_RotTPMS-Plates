%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Title: Isogeometric analysis for Thermo-mechanical buckling of RotTPMS plates with Quasi-3D six-variable plate model %%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% ! Please reference to paper: ............................................
% ! This work can be used, modified, and shared under the MIT License
% ! This work can be found in https://github.com/SIMOGroup/Ther-Mech-Buck-RotTPMS-Plates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% =========================== Initialization =============================
tic
addpath(genpath('./'));

% clc
clear all
% close all
format long

%% ============================ Plate geometry ============================
% === Physical geometric properties ===
Plate.geo.L = 1;
Plate.geo.W = Plate.geo.L;
Plate.geo.h = Plate.geo.L/(10);
% Plate.geo.h = Plate.geo.W/20;

% === Plate theory ===
% [1]: fz = 0, [2]: Reddy (1984), [3]: Shimpi (2002), [4]: Nguyen-Xuan (2013), [5]  Nguyen (2017)
% [6]: Nguyen (2016), [7]: Thai (2014), [8]: Touratier (1991), [9]: Hoang (2023), [10]: Tran (2025)
Plate.theory.shear_func = 2;
% [1]: f'z, [2]: 3/20*f'z, [3]: 1/8*f'z, [4]: 1/12*f'z, [5]: 1/20*f'z, [6]: f'z + 1, [7]: 1 - f'z
Plate.theory.stretch_func = 4;

% === Define NURBS functions ===
IGA.NURBS.deg = 3; % Degree of basis functions
IGA.NURBS.ref = 8; % Number of mesh refinement

%% ============================ Material ==================================
% === Base matrix material ===
% [1]: ZrO2, [2]: Al2O3, [3]: Si3N4, [4]: Ti-6Al-4V, [5]: SUS304, [6]: Ni
Plate.mat.type = 5;
Plate.mat.T_0 = 300;

% === Porous material ===
% [1]: Primitive, [2]: Gyroid, [3]: IWP
Plate.por_mat.type = 1;
Plate.por_mat.RD = 0.35;
% Rxyz order with Extrinsic rotation, Counter-clockwise angles
Plate.por_mat.alpha = [0, 0, 0];
% Plate.por_mat.alpha = [0, 0, 45];
% Plate.por_mat.alpha = [90, -45, 90];
% Plate.por_mat.alpha = [90, -45, 135];
% Plate.por_mat.alpha = [45, -35.264, 0];

load_dir = 2; alpha_AB = 90; theta_z = 45; 
if load_dir == 1
    dir_vec = [0; tan(deg2rad(alpha_AB)); 1];  % Curve A: wrt OYZ plane
%     alpha_x = alpha_AB; alpha_y = 0;
elseif load_dir == 2
    dir_vec = [tan(deg2rad(alpha_AB)); 1; 1];  % Curve B: wrt OXC plane where C is the center point of 1/4 circle in OYZ plane
%     alpha_x = 45; alpha_y = rad2deg(- atan2(tan(deg2rad(alpha_AB)),sqrt(2)));
end

[alpha_x, alpha_y] = compute_Rotation_Angles_xy(dir_vec);
alpha_z = compute_Rotation_Angle_z(alpha_x, alpha_y, 3);  % [1]: X' _|_ OYZ, [2]: [1;1;1] _|_ OYZ, [3]: Depends on subspaces
Plate.por_mat.alpha = [alpha_x, alpha_y, alpha_z + theta_z];
            
% === Temperature distribution ===
% [1]: Uniform (T = T1), [2]: Linear
Plate.temp_dis.type = 1;
Plate.temp_dis.T1 = 300 + 1e-4;  % Bot
Plate.temp_dis.T2 = 300;  % Top

%% ========================== Problem type ================================
% [2]: Vibration, [3]: Thermal Buckling, [4]: Thermo-Mechanical buckling
Plate.prob.type = 3;
switch Plate.prob.type
    case 2  % Thermal free vibration
        IGA.result.norm_method = 1;
        IGA.result.nmode = 1;
    case 3  % Thermal Buckling
        IGA.result.norm_method = 1;
        IGA.result.nmode = 1;
    case 4  % Thermo-Mechanical buckling
        Plate.prob.buck_case = 1;  % Buckling load case [1]: X, [2]: Y, [3]: X+Y
        IGA.result.norm_method = 1;
        IGA.result.nmode = 1;
end

%% ========================= Boundary condition ===========================
% [Left, Lower, Right, Upper]
% [1]: Clamped (C), [2]: Hinged (freely movable) supported (S), [3]: Free (F), [4]: Soft (quasi-movable) supported (S2), [5]: Hard (immovable) supported (S3)
% Plate.bc.bc_case = [1,1,1,1];
Plate.bc.bc_case = [4,4,4,4];

%% ======================== Material matrices =============================
[Plate.mat_mat.D, Plate.mat_mat.I, Plate.mat_mat.S_th] = cal_Material_Matrices_2D_6dof_RotTPMS(Plate);

%% =============================== IGA mesh ===============================
% === Generate NURBS mesh and properties ===    
IGA.NURBS = gen_rec_surf(Plate, IGA.NURBS);
IGA.NURBS = gen_Ien_Inn_surf(IGA.NURBS);
IGA.NURBS = gen_Iee_Ine_surf(IGA.NURBS);
IGA.NURBS = gen_FE_approx_surf(IGA.NURBS);

% Options: all, physical_shape, physical_surf, physical_element, control_net, plot_physical_element, note_physical_element, note_physical_edge, plot_control_net, note_control_net
% plot_NURBS_surf(IGA.NURBS, 'physical_surf')
% plot_NURBS_surf(IGA.NURBS, 'physical_shape', 'control_net')
% plot_NURBS_surf(IGA.NURBS, 'plot_physical_domain', 'control_net')

% === NURBS properties ===
IGA.NURBS.nsd   = 2;                                                             % Number of spatial dimension
IGA.NURBS.nnode = IGA.NURBS.mcp * IGA.NURBS.ncp;                                 % Number of control point
IGA.NURBS.nshl  = (IGA.NURBS.p + 1) * (IGA.NURBS.q + 1);                         % Number of local shape functions (= degree + 1 per element due to k refinement)
IGA.NURBS.nel   = (IGA.NURBS.mcp - IGA.NURBS.p) * (IGA.NURBS.ncp - IGA.NURBS.q); % Number of element

% === IGA properties ===
IGA.params.ndof   = 6;                                                           % Number of dofs of a control point
IGA.params.sdof   = IGA.NURBS.nnode * IGA.params.ndof;                           % Total number of dofs of the structure
IGA.params.nGauss = IGA.NURBS.p + 1;                                             % Number of gauss point in integration

%% ========================= IGA for linear geometry ======================
% === Building global matrices === 
IGA.result.K = cal_Stiffness_Matrices_2D_6dof(IGA,Plate);                    % Stiffness
IGA.result.K_th = cal_Thermal_Stiffness_Matrices_2D_6dof(IGA,Plate);         % Thermal Stiffness
switch Plate.prob.type
    case 2  % Thermal free vibration
        IGA.result.M = cal_Mass_Matrices_2D_6dof(IGA,Plate);                 % Mass
    case 4  % Thermo-mechanical buckling
        IGA.result.Kg = cal_Geometric_Stiffness_Matrices_2D_6dof(IGA,Plate); % Geometric Stiffness     
end

% === Imposing boundary conditions ===
[IGA.params.bcdof, IGA.params.bcval] = cal_bcdof_2D_6dof(IGA,Plate);
IGA.params.fdof = setdiff((1:IGA.params.sdof)', IGA.params.bcdof');  % Free dofs

%% === Solving weak form ===
switch Plate.prob.type
    case 2  % Thermal free vibration
        disp("--- Thermal free vibration responses ---")
        
        % --- Calculations ---
        sdof = IGA.params.sdof; fdof = IGA.params.fdof;
        KK = full(IGA.result.K); KK_th = full(IGA.result.K_th); MM = full(IGA.result.M);
        [ModeShape, Lambda] = eig(KK(fdof, fdof) + KK_th(fdof, fdof), MM(fdof, fdof));
        
        [IGA.result.Lambda, sort_index] = sort(diag(Lambda),'ascend');
        IGA.result.ModeShape = zeros(sdof, length(fdof));
        IGA.result.ModeShape(fdof,:) = ModeShape(:, sort_index);
        IGA.result.Lambda = IGA.result.Lambda(IGA.result.Lambda > 1e-4 & IGA.result.Lambda < inf);
        IGA.result.ModeShape = IGA.result.ModeShape(:, IGA.result.Lambda > 1e-4 & IGA.result.Lambda < inf);
        
        % --- Vibration frequencies ---
        nmode = IGA.result.nmode;
        L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h;
        [E_m, nu_m, rho_m, ~, ~] = compute_temperature_dependent_material(Plate.mat.type, 300);
%         [E_m, nu_m, rho_m, ~] = compute_basic_material(Plate.mat.type);
        I_m = h*rho_m; D_m = E_m*h^(3)/(12*(1-nu_m^2));
        switch IGA.result.norm_method
            case 1
                Lambda_norm = (IGA.result.Lambda(1:nmode)*rho_m*L^4*h/D_m).^0.25;
            case 2
                norm = h*sqrt(rho_m/E_m);
                Lambda_norm = sqrt(IGA.result.Lambda(1:nmode))*norm;
            case 3
                norm = L*sqrt(rho_m*(1-nu_m^(2))/E_m);
                Lambda_norm = sqrt(IGA.result.Lambda(1:nmode))*norm;
            case 4
                norm = L^2/h*sqrt(rho_m/E_m);
                Lambda_norm = sqrt(IGA.result.Lambda(1:nmode))*norm;
            case 5
                norm = W^(2)/pi^(2)*(I_m/D_m)^(1/2);
                Lambda_norm = sqrt(IGA.result.Lambda(1:nmode))*norm;
            case 6
                norm = W^(2)/h*(rho_m*(1-nu_m^2)/E_m)^(1/2);
                Lambda_norm = sqrt(IGA.result.Lambda(1:nmode))*norm;
            case 7
                norm = 1/2/pi;
                Lambda_norm = sqrt(IGA.result.Lambda(1:nmode))*norm;
            case 8 
                norm = L^2/h*sqrt(rho_m/E_m);
                Lambda_norm = sqrt(IGA.result.Lambda(1:nmode))*norm;
        end
        disp("> Normalized frequencies = " + sprintf('%.4f, ', Lambda_norm))
        
%         % --- Modeshape displacement NURBS mesh ---
%         for i_mode = 1:nmode
%             IGA_mode = IGA; Plate_mode = Plate;
%             IGA_mode.result.U = IGA.result.ModeShape(:, i_mode); Plate_mode.prob.q_uniform = 1; IGA_mode.result.norm_method = 0;
%             NURBS_disp_norm = cal_displacement_2D_6dof(IGA_mode,Plate_mode);
%             NURBS_disp = IGA.NURBS; plot_size_factor = 1/max(abs(NURBS_disp_norm(:)));
%             NURBS_disp.CP(:,:,3) = NURBS_disp.CP(:,:,3) + plot_size_factor * NURBS_disp_norm;
%             NURBS_disp = gen_Iee_Ine_surf(NURBS_disp); NURBS_disp = gen_FE_approx_surf(NURBS_disp);
%             plot_NURBS_surf(NURBS_disp, 'physical_surf', 'plot_physical_element', 'plot_control_net')
%         end
        
        clear sdof fdof nmode KK KK_th MM ModeShape Lambda sort_index
        clear L W h E_m nu_m rho_m I_m D_m 
        clear i_mode IGA_mode Plate_mode plot_size_factor
            case 3  % Thermal buckling
        disp("--- Thermal buckling responses ---")
        
        % --- Calculations --- 
        sdof = IGA.params.sdof; fdof = IGA.params.fdof;
        KK = full(IGA.result.K); KK_th = full(IGA.result.K_th);
        [ModeShape, Lambda] = eig(KK(fdof, fdof), -KK_th(fdof, fdof));
        
        [IGA.result.Lambda, sort_index] = sort(diag(Lambda),'ascend');
        IGA.result.ModeShape = zeros(sdof, length(fdof));
        IGA.result.ModeShape(fdof,:) = ModeShape(:, sort_index);
        IGA.result.Lambda = IGA.result.Lambda(IGA.result.Lambda > 1e-4 & IGA.result.Lambda < inf);
        IGA.result.ModeShape = IGA.result.ModeShape(:, IGA.result.Lambda > 1e-4 & IGA.result.Lambda < inf);
        
        % --- Critical buckling parameters ---
        nmode = IGA.result.nmode;
        L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h;
        [E_m, nu_m, rho_m, alpha_m, ~] = compute_temperature_dependent_material(Plate.mat.type, 300);
%         [E_m, nu_m, rho_m, alpha_m] = compute_basic_material(Plate.mat.type);
        I_m = h*rho_m; D_m = E_m*h^(3)/(12*(1-nu_m^2));
        switch IGA.result.norm_method
            case 1
                norm = (Plate.temp_dis.T1 - Plate.mat.T_0);
                Lambda_norm = IGA.result.Lambda(1:nmode)*norm;
            case 2
                norm = alpha_m * (Plate.temp_dis.T1 - Plate.mat.T_0);
                Lambda_norm = IGA.result.Lambda(1:nmode)*norm;
        end
        disp("> Normalized critical buckling temperature rise = " + sprintf('%.4f, ', Lambda_norm))
        
%         % --- Modeshape displacement NURBS mesh ---
%         for i_mode = 1:nmode
%             IGA_mode = IGA; Plate_mode = Plate;
%             IGA_mode.result.U = IGA.result.ModeShape(:, i_mode); Plate_mode.prob.q_uniform = 1; IGA_mode.result.norm_method = 0;
%             NURBS_disp_norm = cal_displacement_2D_6dof(IGA_mode,Plate_mode);
%             NURBS_disp = IGA.NURBS; plot_size_factor = 1/max(abs(NURBS_disp_norm(:)));
%             NURBS_disp.CP(:,:,3) = NURBS_disp.CP(:,:,3) + plot_size_factor * NURBS_disp_norm;
%             NURBS_disp = gen_Iee_Ine_surf(NURBS_disp); NURBS_disp = gen_FE_approx_surf(NURBS_disp);
%             plot_NURBS_surf(NURBS_disp, 'physical_surf', 'plot_physical_element', 'plot_control_net')
%         end
        
        clear sdof fdof nmode KK KK_th ModeShape Lambda sort_index
        clear L W h E_m nu_m rho_m I_m D_m 
        clear i_mode IGA_mode Plate_mode plot_size_factor
        
    case 4  % Thermo-mechanical buckling
        disp("--- Thermo-mechanical buckling responses ---")
        
        % --- Calculations --- 
        sdof = IGA.params.sdof; fdof = IGA.params.fdof;
        KK = full(IGA.result.K); KK_th = full(IGA.result.K_th); KKg = full(IGA.result.Kg);
        [ModeShape, Lambda] = eig(KK(fdof, fdof)+KK_th(fdof, fdof), KKg(fdof, fdof));
        
        [IGA.result.Lambda, sort_index] = sort(diag(Lambda),'ascend');
        IGA.result.ModeShape = zeros(sdof, length(fdof));
        IGA.result.ModeShape(fdof,:) = ModeShape(:, sort_index);
        IGA.result.Lambda = IGA.result.Lambda(IGA.result.Lambda > 1e-4 & IGA.result.Lambda < inf);
        IGA.result.ModeShape = IGA.result.ModeShape(:, IGA.result.Lambda > 1e-4 & IGA.result.Lambda < inf);
        
        % --- Critical buckling parameters ---
        nmode = IGA.result.nmode;
        L = Plate.geo.L; W = Plate.geo.W; h = Plate.geo.h;
        [E_m, nu_m, rho_m, alpha_m, ~] = compute_temperature_dependent_material(Plate.mat.type, 300);
%         [E_m, nu_m, rho_m, ~] = compute_basic_material(Plate.mat.type);
        I_m = h*rho_m; D_m = E_m*h^(3)/(12*(1-nu_m^2));
        switch IGA.result.norm_method
            case 1
                norm = W^2/(pi^2*D_m);
                Lambda_norm = IGA.result.Lambda(1:nmode)*norm;
            case 2
                norm = L^2/D_m;
                Lambda_norm = IGA.result.Lambda(1:nmode)*norm;
        end
        disp("> Normalized critical buckling load = " + sprintf('%.4f, ', Lambda_norm))
        
%         % --- Modeshape displacement NURBS mesh ---
%         for i_mode = 1:nmode
%             IGA_mode = IGA; Plate_mode = Plate;
%             IGA_mode.result.U = IGA.result.ModeShape(:, i_mode); Plate_mode.prob.q_uniform = 1; IGA_mode.result.norm_method = 0;
%             NURBS_disp_norm = cal_displacement_2D_6dof(IGA_mode,Plate_mode);
%             NURBS_disp = IGA.NURBS; plot_size_factor = 1/max(abs(NURBS_disp_norm(:)));
%             NURBS_disp.CP(:,:,3) = NURBS_disp.CP(:,:,3) + plot_size_factor * NURBS_disp_norm;
%             NURBS_disp = gen_Iee_Ine_surf(NURBS_disp); NURBS_disp = gen_FE_approx_surf(NURBS_disp);
%             plot_NURBS_surf(NURBS_disp, 'physical_surf', 'plot_physical_element', 'plot_control_net')
%         end
        
        clear sdof fdof nmode KK KK_th KKg ModeShape Lambda sort_index
        clear L W h E_m nu_m rho_m I_m D_m 
        clear i_mode IGA_mode Plate_mode plot_size_factor
end
toc
