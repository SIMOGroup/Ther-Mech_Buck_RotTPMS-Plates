function [E, nu, rho, alpha] = compute_basic_material(mat_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Material library %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
% [E (Pa), nu, rho (kg/m^3), alpha (1/K)]
% [1]: Steel, [2]: Alumium, [3]: Titanium, [4]: Copper, [5]: Brass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Unit: m - N - kg - s^2 - Pa
switch mat_type
    case 1  % Steel
        E = 200e9; nu = 0.3; rho = 8000; alpha = 1e-6;
    case 2  % Alumium
        E = 70e9; nu = 0.3; rho = 2702; alpha = 1e-6;
    case 3  % Titanium
        E = 106e9; nu = 0.3; rho = 4510; alpha = 1e-6;
    case 4  % Copper
        E = 117e9; nu = 0.3; rho = 8960; alpha = 1e-6;
    case 5  % Brass
        E = 90e9; nu = 0.3; rho = 8730; alpha = 1e-6;
    case 6  % Test
        E = 1e6; nu = 0.3; rho = 8000; alpha = 1e-6;
end
end
