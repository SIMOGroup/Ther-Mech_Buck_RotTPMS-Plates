function [bcdof, bcval] = cal_bcdof_2D_6dof(IGA,Plate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate dofs of boundary conditions with Quasi-3D 6-variable plate theory %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from IGA
mcp = IGA.NURBS.mcp; ncp = IGA.NURBS.ncp;
ndof = IGA.params.ndof;

%% Used parameters from Plate
bc_case = Plate.bc.bc_case;

%% ===== Initial boundary condition arrays =====
bcdof = []; bcval = [];

%% ===== Control points on boundaries =====
CP_peri_1 = {[ 1: mcp: mcp*(ncp-1)+1];  % Left   % Perimeter
             [                1: mcp];  % Lower
             [     mcp: mcp: mcp*ncp];  % Right
             [mcp*(ncp-1)+1: mcp*ncp]}; % Upper
       
CP_peri_2 = {[     2: mcp: mcp*(ncp-1)+2];  % Left   % 01 control point inside
             [              mcp+1: 2*mcp];  % Lower
             [     mcp-1: mcp: mcp*ncp-1];  % Right
             [mcp*(ncp-2)+1: mcp*(ncp-1)]}; % Upper

%% ===== Dofs of control points on boundaries (dofs = u_{0}, v_{0}, w_{0}, \beta_{x}, \beta_{y}, \beta_{z}) =====
for i_bc = 1:length(bc_case)
    switch bc_case(i_bc)
        case 1  % Clamped (C)
            for idof = [1, 2, 3, 4, 5, 6]  % Constraint u_{0}, v_{0}, w_{0}, \beta_{x}, \beta_{y}, \beta_{z} at perimeter
                bcdof = [bcdof ndof*(CP_peri_1{i_bc}-1)+idof];
            end

            for idof = [3, 6]  % Constraint dw_{0}dx, dw_{0}dy, d\beta_{z}dx, d\beta_{z}dy at perimeter (with a fine mesh)
                bcdof = [bcdof ndof*(CP_peri_2{i_bc}-1)+idof];
            end
        case 2  % Hinged (freely movable) supported
            for idof = [3, 6]  % Constraint w_{0} and \beta_{z} at perimeter
                bcdof = [bcdof ndof*(CP_peri_1{i_bc}-1)+idof];
            end
            
        case 4  % Soft (quasi-movable) supported (S2)
            for idof = [3, 6]  % Constraint w_{0} and \beta_{z} at perimeter
                bcdof = [bcdof ndof*(CP_peri_1{i_bc}-1)+idof];
            end
            
            if i_bc == 2 || i_bc == 4
                for idof = [1, 4]  % Constraint u_{0}, \beta_{x} at perimeter on Upper and Lower
                    bcdof = [bcdof ndof*(CP_peri_1{i_bc}-1)+idof];
                end
            elseif i_bc == 1 || i_bc == 3
                for idof = [2, 5]  % Constraint v_{0}, \beta_{y} at perimeter on Left and Right
                    bcdof = [bcdof ndof*(CP_peri_1{i_bc}-1)+idof];
                end
            end
        case 5  % Hard (immovable) supported (S3)
            for idof = [3, 6]  % Constraint w_{0} and \beta_{z} at perimeter
                bcdof = [bcdof ndof*(CP_peri_1{i_bc}-1)+idof];
            end
            
            for idof = [1, 2, 4, 5]  % Constraint u_{0}, v_{0}, \beta_{x}, \beta_{y} at all perimeters
                bcdof = [bcdof ndof*(CP_peri_1{i_bc}-1)+idof];
            end
    end
end
bcval = [bcval, zeros(1,length(bcdof))];

end
