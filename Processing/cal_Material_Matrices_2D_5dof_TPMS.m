function [Db, Ds, I, S_th] = cal_Material_Matrices_2D_5dof_TPMS(Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate material matrices for TPMS with five-variable plate theory %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from Plate
h = Plate.geo.h;
shear_func = Plate.theory.shear_func;
mat_type = Plate.mat.type; T0 = Plate.mat.T_0;
porous_type = Plate.por_mat.type; RD = Plate.por_mat.RD;
temp_dis = Plate.temp_dis.type; T1 = Plate.temp_dis.T1 ; T2 = Plate.temp_dis.T2;

%% ===== Initial material matrices =====
for i = 3:-1:1
    for j = 3:-1:1
        Db{i,j} = zeros(3,3);
    end
end
Ds = zeros(2,2);
for i = 3:-1:1
    for j = 3:-1:1
        I{i,j} = zeros(3,3);
    end
end

for i = 3:-1:1
    for j = 3:-1:1
        S_th{i,j} = zeros(2,2);
    end
end

%% ===== Gauss integration =====
[Pg, Wg] = gen_Gauss_point(1,-1,20) ;

%% ===== Define material properties =====
% --- TPMS parameters ---
switch porous_type
    case 1 % Primitive
        k_e = 0.25; C_1e = 0.317; n_1e = 1.264; n_2e = 2.006;
        k_g = 0.25; C_1g = 0.705; n_1g = 1.189; n_2g = 1.715;
        k_nu = 0.55; a_1 = 0.314; b_1 = -1.004; a_2 = 0.152;
    case 2 % Gyroid
        k_e = 0.45; C_1e = 0.596; n_1e = 1.467; n_2e = 2.351;
        k_g = 0.45; C_1g = 0.777; n_1g = 1.544; n_2g = 1.982;
        k_nu = 0.50; a_1 = 0.192; b_1 = -1.349; a_2 = 0.402;
    case 3 % IWP
        k_e = 0.35; C_1e = 0.597; n_1e = 1.225; n_2e = 1.782;
        k_g = 0.35; C_1g = 0.529; n_1g = 1.287; n_2g = 2.188;
        k_nu = 0.13; a_1 = 2.597; b_1 = -0.157; a_2 = 0.201;
end
C_2e = (C_1e*k_e^(n_1e) - 1)/(k_e^(n_2e) - 1); C_3e = 1 - C_2e;
C_2g = (C_1g*k_g^(n_1g) - 1)/(k_g^(n_2g) - 1); C_3g = 1 - C_2g;
d_1 = 0.3 - a_1*exp(b_1*k_nu); b_2 = - a_2*(k_nu + 1); d_2 = 0.3 - a_2*(1)^2 - b_2(1);

%% ===== Calculate material matrices =====
z_lower = -h/2; z_upper = h/2;

for iGauss = 1:size(Wg,1)  % Loop over the integration points
    % --- Gauss points & weights ---
    gpt = Pg(iGauss);
    gwt = Wg(iGauss);
    
    % --- Map the point to global ---
    gpt = (z_upper-z_lower)/2*gpt + (z_upper+z_lower)/2;
    gwt = (z_upper-z_lower)/2*gwt;
	
    % --- Temperature distribution ---
    switch temp_dis
        case 1
            T = T1;
        case 2
            T = T1 + (1/2 + gpt/h)*(T2 - T1);
    end
    Delta_T = T - T0;
    
    % --- Material properties ---
    [E_s, nu_s, rho_s, alpha_s, ~] = compute_temperature_dependent_material(mat_type, T);
%     [E_s, nu_s, rho_s, alpha_s] = compute_basic_material(mat_type);
    G_s = E_s / (2*(1+nu_s));
    
    % --- Porous material properties ---
    e = (RD <= k_e) * (C_1e*RD^(n_1e)) + ...
        (RD >  k_e) * (C_2e*RD.^(n_2e) + C_3e);
    g = (RD <= k_g) * (C_1g*RD^(n_1g)) + ...
        (RD >  k_g) * (C_2g*RD.^(n_2g) + C_3g);
    nu =(RD <= k_nu) * (a_1.*exp(b_1*RD) + d_1) + ...
        (RD >  k_nu) * (a_2*RD^(2) + b_2*RD + d_2);
    E = E_s*e;
    G = G_s*g;
    rho = rho_s * RD;
    alpha = alpha_s;
    
    % --- Consitutive matrix ---
    C_11 = E / (1-nu^(2)); C_12 = (E*nu) / (1-nu^(2)); C_44 = G;
    Qb = [C_11, C_12,    0; ...
          C_12, C_11,    0; ...
             0,    0, C_44];
    Qs = [C_44,    0; ...
             0, C_44];
    
    % --- Shear deformation function ---
    [fz, dfz, ~] = compute_shear_deformation_function(gpt,h,shear_func);
    
    % --- Bending stiffness matrix ---
    h_e = [1 gpt fz];
    for i = 1:3
        for j = 1:3
            Db{i,j} = Db{i,j} + h_e(i)*h_e(j)*Qb*gwt;
        end
    end
    
    % --- Shear stiffness matrix ---
    Ds = Ds + dfz^(2)*Qs*gwt;
    
    % --- Inertia matrix ---
    h_u = [1 gpt fz];
    for i = 1:3
        for j = 1:3
            I{i,j} = I{i,j} + h_u(i)*h_u(j)*rho*eye(3)*gwt;
        end
    end
    
    % --- In-plane thermal stress matrix --- 
    S_th_0 = - alpha * Delta_T * [C_11 + C_12,           0; ...
                                            0, C_12 + C_11];
    for i = 1:3
        for j = 1:3
            S_th{i,j} = S_th{i,j} + h_u(i)*h_u(j) * S_th_0 * gwt;
        end
    end
    
    % --- Neglecting in-plane displacement in buckling mode ---
%     S_th{1,2}(:,:) = 0; S_th{1,3}(:,:) = 0; 
%     S_th{2,1}(:,:) = 0; S_th{2,2}(:,:) = 0; S_th{2,3}(:,:) = 0; 
%     S_th{3,1}(:,:) = 0; S_th{3,2}(:,:) = 0; S_th{3,3}(:,:) = 0;
end
end
