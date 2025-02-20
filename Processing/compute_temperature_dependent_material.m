function [E, nu, rho, alpha, k] = compute_temperature_dependent_material(mat_type, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Temperature-dependent material library %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
% [E (Pa), nu, rho (kg/m^3), alpha (1/K), k (W/mK)]
% [1]: ZrO2, [2]: Al2O3, [3]: Si3N4, [4]: Ti-6Al-4V, [5]: SUS304, [6]: Ni, [7]: PmPV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Unit: m - N - kg - s^2 - Pa
switch mat_type  % (P_0, P_-1, P_1, P_2, P_3)
    case 1  % ZrO2 (Zirconia)
        P_E = [244.27e9, 0, -1.371e-3,  1.214e-6, -3.681e-10];
        P_nu = [0.2882, 0,  1.133e-4,        0, 0];
        P_rho = [3000, 0, 0, 0, 0];
        P_alpha = [12.766e-6, 0, -1.491e-3,  1.006e-5, -6.778e-11];       
        P_k = [   1.700,       0,  1.276e-4, 6.648e-8,          0];

    case 2  % Al2O3 (Aluminium oxide)
        P_E = [349.55e9, 0, -3.853e-4,  4.027e-7, -1.673e-10];
        P_nu = [0.2600, 0,         0,        0, 0];
        P_rho = [3800, 0, 0, 0, 0];
        P_alpha = [ 6.827e-6, 0,  1.838e-4,         0,          0];
        P_k = [ -14.087, -1123.6, -6.227e-3,        0,          0];
        
    case 3  % Si3N4 (Silicon nitride)
        P_E = [348.43e9, 0, -3.070e-4,  2.160e-7, -8.946e-11];
        P_nu = [0.2400, 0,         0,        0, 0];
        P_rho = [2370, 0, 0, 0, 0];
        P_alpha = [ 5.872e-6, 0,  9.095e-4,         0,          0];
        P_k = [  13.723,       0, -1.032e-3, 5.466e-7, -7.876e-11];
        
    case 4  % Ti-6Al-4V
        P_E = [122.56e9, 0, -4.586e-4,         0,          0];
        P_nu = [0.2884, 0,  1.121e-4,        0, 0];
        P_rho = [4420, 0, 0, 0, 0];
        P_alpha = [ 7.579e-6, 0,  6.638e-4, -3.147e-6,          0];
        P_k = [   1.000,       0,  1.704e-2,        0,          0];
        
    case 5  % SUS304 (Stainless steel)
        P_E = [201.04e9, 0,  3.079e-4, -6.534e-7,          0];
        P_nu = [0.3262, 0, -2.002e-4, 3.797e-7, 0];
        P_rho = [8166, 0, 0, 0, 0];
        P_alpha = [12.330e-6, 0,  8.086e-4,         0,          0];
        P_k = [  15.379,       0, -1.264e-3, 2.092e-6, -7.223e-10];
        
    case 6  % Ni (Nickel)
        P_E = [223.95e9, 0, -2.794e-4,  3.998e-9,          0];
        P_nu = [0.3100, 0,         0,        0, 0];
        P_rho = [3000, 0, 0, 0, 0];
        P_alpha = [ 9.921e-6, 0,  8.705e-4,         0,          0];
        if T <= 635
            P_k = [ 187.660,       0, -2.869e-3, 4.005e-6, -1.983e-9 ];
        else
            P_k = [  58.754,       0, -4.614e-4, 6.657e-7, -1.523e-10];
        end
        
    case 7  % Poly[(m-phenylenevinylene)-co-(2,5-dioctoxy-p-phenylenevinylene)] (PmPV)
        T_0 = 300;
        E = (2.1 - 0.0047*(T - T_0))*1e9; nu = 0.34; rho = 1190;
        alpha = 45*(1 + 0.0005*(T - T_0))*1e-6; k = 0;
        
    case 8  % Steel        
        E = 200e9; nu = 0.3; rho = 8000;
        alpha = 0; k = 0;
        
    case 9  % Test
        E = 1e6; nu = 0.3; rho = 1000;
        alpha = 1e-6; k = 0;
end

if ~exist('E', 'var')
    E = P_E(1) * (P_E(2)*1/T + 1 + P_E(3)*T + P_E(4)*T^2 + P_E(5)*T^3);
    nu = P_nu(1) * (P_nu(2)*1/T + 1 + P_nu(3)*T + P_nu(4)*T^2 + P_nu(5)*T^3);
    rho = P_rho(1) * (P_rho(2)*1/T + 1 + P_rho(3)*T + P_rho(4)*T^2 + P_rho(5)*T^3);
    alpha = P_alpha(1) * (P_alpha(2)*1/T + 1 + P_alpha(3)*T + P_alpha(4)*T^2 + P_alpha(5)*T^3);
    k = P_k(1) * (P_k(2)*1/T + 1 + P_k(3)*T + P_k(4)*T^2 + P_k(5)*T^3);
end

end
