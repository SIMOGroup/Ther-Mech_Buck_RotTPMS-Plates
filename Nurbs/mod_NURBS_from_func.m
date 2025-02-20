function NURBS = mod_NURBS_from_func(NURBS, z_func, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate NURBS surf for 2D rectangle with function of z-coordinate %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from NURBS
p = NURBS.p; uKnot = NURBS.uKnot; mcp = NURBS.mcp;
q = NURBS.q; vKnot = NURBS.vKnot; ncp = NURBS.ncp;
CP = NURBS.CP;

%% ===== Modify NURBS =====
% --- Define approximation shape matrix A coordinates ---
npoint = 101;
FE_grid = zeros(npoint, npoint, 3);

A = zeros(numel(FE_grid(:,:,1)), numel(CP(:,:,1)));
for i_point_y = 1:npoint
    v = vKnot(q+1) + (vKnot(ncp+1) - vKnot(q+1))/(npoint - 1) * (i_point_y-1);
    nj = find_knot_span(v,vKnot,ncp);
    for i_point_x = 1:npoint
        u = uKnot(p+1) + (uKnot(mcp+1) - uKnot(p+1))/(npoint - 1) * (i_point_x-1);
        ni = find_knot_span(u,uKnot,mcp);

        idx_point = (i_point_y-1)*npoint + i_point_x;
        Nu = eval_basis_func(ni,p,u,uKnot); 
        Nv = eval_basis_func(nj,q,v,vKnot);

        SumNw = 0;
        for j = 0:q
            for i = 0:p
                SumNw = Nu(i+1)*Nv(j+1)*CP(ni-p+i,nj-q+j,4) + SumNw;
            end
        end

        for j = 0:q
            for i = 0:p
                idx_ctrl = (nj-q+j-1)*(mcp) + ni-p+i;
                A(idx_point, idx_ctrl) = Nu(i+1)*Nv(j+1)*CP(ni-p+i,nj-q+j,4)/SumNw;
            end
        end
    end
end

for i_dim = 1:3
    FE_grid(:,:,i_dim) = reshape(A*reshape(CP(:,:,i_dim),[],1),npoint,npoint);
end

% --- Modify z coordinate ---
P_true = reshape(z_func(FE_grid(:,:,1),FE_grid(:,:,2)),[],1);
CP_approx = pinv(A,1e-8)*P_true;
P_approx = A*CP_approx;

P_true = reshape(P_true,npoint,npoint);
CP_approx = reshape(CP_approx,mcp,ncp);
P_approx = reshape(P_approx,npoint,npoint);

NURBS.CP(:,:,3) = CP_approx;

% --- Metrics ---
NURBS.Approx.RMSE = sqrt(mean((P_true(:) - P_approx(:)).^2));
if ~isempty(find(P_true(:)==0, 1))
    NURBS.Approx.MAPE = nan;
else
    NURBS.Approx.MAPE = mean(abs((P_true(:) - P_approx(:))./(P_true(:))))*100;
end

if mean(P_true(:) - mean(P_true(:))) == 0
    NURBS.Approx.R2 = 1;
else
    NURBS.Approx.R2 = 1 - mean((P_true(:) - P_approx(:)).^2) / mean((P_true(:) - mean(P_true(:))).^2);
end

disp("--- Approximation metrics ---")
disp("> RMSE = " + sprintf('%.4e', NURBS.Approx.RMSE))
disp("> MAPE = " + sprintf('%.4f', NURBS.Approx.MAPE) + "%")
disp("> R2 = " + sprintf('%.4f', NURBS.Approx.R2))

%% ===== NURBS properties =====
NURBS = gen_Ien_Inn_surf(NURBS);
NURBS = gen_Iee_Ine_surf(NURBS);
NURBS = gen_FE_approx_surf(NURBS);

NURBS.nsd   = 2;                                               % Number of spatial dimension
NURBS.nnode = NURBS.mcp * NURBS.ncp;                           % Number of control point
NURBS.nshl  = (NURBS.p + 1) * (NURBS.q + 1);                   % Number of local shape functions (= degree + 1 per element due to k refinement)
NURBS.nel   = (NURBS.mcp - NURBS.p) * (NURBS.ncp - NURBS.q);   % Number of element

% --- Calculate RD_avg ---
uKnot = NURBS.uKnot; vKnot = NURBS.vKnot;
Inn = NURBS.Inn; Ien = NURBS.Ien; 
nGauss = 5; [Pg, Wg] = gen_Gauss_point(1,-1,nGauss);

RD_I = NURBS.CP(:,:,3); RD_total = 0; Area = 0;
for iel = 1:NURBS.nel  % Loop over the elements
    % --- Element parameters ---
    sctr = Ien(iel,:);  % Control points indexes
    ni = Inn(Ien(iel,1),1);  % Index of the element in parametric domain
    nj = Inn(Ien(iel,1),2);
    
    % --- Gauss integration ---
    if abs(uKnot(ni)-uKnot(ni+1)) > 1e-8 && abs(vKnot(nj)-vKnot(nj+1)) > 1e-8  % Check if the current element has nonzero area in the parametric domain
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
                [N, ~, ~, ~, ~, detJ1] = Kine_Shape_2nd_2D(NURBS,ni,nj,gpt_xi,gpt_eta);
                
                % Material matrices
                RD_total = RD_total + N'*RD_I(sctr)'*gwt*detJ1;
                Area = Area + gwt*detJ1;
            end
        end
    end
end
NURBS.RD_avg = RD_total / Area;

%% ===== Plot reconstruction surf =====
%% Check plot options
plot_original_surf = false; plot_reconstructed_surf = false; plot_error_surf = false;

for i = 1:length(varargin)
    switch varargin{i}
        case 'all'
            plot_original_surf = true;
            plot_reconstructed_surf = true;
            plot_error_surf = true;
        case 'plot_original_surf'
            plot_original_surf = true;
        case 'plot_reconstructed_surf'
            plot_reconstructed_surf = true;
        case 'plot_error_surf'
            plot_error_surf = true;
        otherwise
            error('Unknown keyword: %s', varargin{i});
    end
end

% --- Original and reconstructed surf ---
if plot_original_surf
    figure('Color','w', 'Units','normalized', 'Outerposition',[0 0 0.6 0.9]);
    set(axes, 'Position', [0.1, 0.1, 0.8, 0.8]);
    hold on

    title('\bf{Original FE Surface}', 'Units','normalized', 'Position',[0.5, 1.05, 0], 'Interpreter','latex', 'fontsize', 14);
    xlabel('x', 'Interpreter','latex', 'fontsize', 14);
    ylabel('y', 'Interpreter','latex', 'fontsize', 14);
    zlabel('z', 'Interpreter','latex', 'fontsize', 14);

    set(gca, 'TickLabelInterpreter','latex', 'TickLength', [0.02 0.02], 'TickDir','in', 'fontsize', 12)
    set(gca, 'box','on', 'GridLineStyle', 'none', 'GridAlpha', 0.2)
    set(gca, 'Layer', 'Top')
    view(3)
    axis equal

    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [22, 29.7]);
    set(gcf, 'PaperOrientation','landscape');
    
    surf(FE_grid(:,:,1), FE_grid(:,:,2), P_true, 'FaceColor','interp', 'FaceAlpha',1, ...
                                                 'EdgeColor', 'k', 'EdgeAlpha', 0.5, 'LineStyle', '-', 'LineWidth', 1, ...
                                                 'Marker', 'none')
end

if plot_reconstructed_surf
    plot_NURBS_surf(NURBS, 'physical_surf', 'plot_control_net')
    title('\bf{Reconstructed NURBS Surf}', 'Units','normalized', 'Position',[0.5, 1.05, 0], 'Interpreter','latex', 'fontsize', 14);
end
    
% --- Errors ---
if plot_error_surf
    figure('Color','w', 'Units','normalized', 'Outerposition',[0 0 0.6 0.7]);
    set(axes, 'Position', [0.1, 0.1, 0.75, 0.8]);
    hold on

    title('\bf{Reconstructed errors}', 'Units','normalized', 'Position',[0.5, 1.05, 0], 'Interpreter','latex', 'fontsize', 14);
    xlabel('x', 'Interpreter','latex', 'fontsize', 14);
    ylabel('y', 'Interpreter','latex', 'fontsize', 14);
    % zlabel('Error', 'Interpreter','latex', 'fontsize', 14);

    set(gca, 'TickLabelInterpreter','latex', 'TickLength', [0.02 0.02], 'TickDir','in', 'fontsize', 12)
    set(gca, 'box','off', 'GridLineStyle', 'none', 'GridAlpha', 0.2)
    set(gca, 'Layer', 'Top')
    view(2)
    axis equal

    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [22, 29.7]);
    set(gcf, 'PaperOrientation','landscape');

    surf(FE_grid(:,:,1), FE_grid(:,:,2), (P_true - P_approx), ...
        'FaceColor','interp', 'FaceAlpha',1, ...
        'EdgeColor', 'none', 'EdgeAlpha', 0.5, 'LineStyle', '-', 'LineWidth', 1, ...
        'Marker', 'none')
    xlim([min(min(FE_grid(:,:,1))), max(max(FE_grid(:,:,1)))])
    ylim([min(min(FE_grid(:,:,2))), max(max(FE_grid(:,:,2)))])

    cbar = colorbar('Position', [0.88 0.16 0.025 0.68], 'TickLength', 0.02, 'TickLabelInterpreter', 'latex');
    colormap(jet);
    cbar.FontSize = 14;
    cbar.Label.String = 'Error';
    cbar.Label.Interpreter = "latex";
    cbar.Label.FontSize = 14;
end

end
