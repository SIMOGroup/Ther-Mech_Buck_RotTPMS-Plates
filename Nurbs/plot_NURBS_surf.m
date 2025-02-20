function plot_NURBS_surf(NURBS, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot NURBS surf with options %%%
% Author: Kim Q. Tran, H. Nguyen-Xuan
% Contact: CIRTech Institude, HUTECH university, Vietnam
% Email: tq.kim@hutech.edu.vn, ngx.hung@hutech.edu.vn
% ! This work can be used, modified, and shared under the MIT License
% Use note:
% 'all' = 'plot_physical_shape', 'plot_physical_surf', 'plot_physical_element', 'note_physical_element', 'note_physical_edge', 'plot_control_net', 'note_control_net'
% 'physical_shape' = 'plot_physical_shape'
% 'physical_surf' = 'plot_physical_shape', 'plot_physical_surf'
% 'physical_element' = 'plot_physical_shape', 'plot_physical_element', 'note_physical_edge', 'note_physical_element' 
% 'control_net' = 'plot_control_net', 'note_control_net'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Used parameters from NURBS
CP = NURBS.CP; FE_grid = NURBS.FE_grid;
Ine = NURBS.Ine; Iee = NURBS.Iee;

%% Check plot options
plot_physical_shape = false; plot_physical_surf = false;
plot_physical_element = false; note_physical_element = false; note_physical_edge = false; 
plot_control_net = false; note_control_net = false;

for i = 1:length(varargin)
    switch lower(varargin{i})
        case 'all'
            plot_physical_shape = true;
            plot_physical_surf = true;
            plot_physical_element = true;
            note_physical_element = true;
            note_physical_edge = true;            
            plot_control_net = true;
            note_control_net = true;
            
        case 'physical_shape'
            plot_physical_shape = true;
        case 'physical_surf'
            plot_physical_shape = true;
            plot_physical_surf = true;
        case 'physical_element'
            plot_physical_shape = true;
            plot_physical_element = true;
            note_physical_edge = true;
            note_physical_element = true;
        case 'control_net'
            plot_control_net = true;
            note_control_net = true;
            
        case 'plot_physical_surf'
            plot_physical_surf = true;
        case 'plot_physical_element'
            plot_physical_element = true;
        case 'note_physical_element'
            note_physical_element = true;
        case 'note_physical_edge'
            note_physical_edge = true;
        case 'plot_control_net'
            plot_control_net = true;
        case 'note_control_net'
            note_control_net = true;
        otherwise
            error('Unknown keyword: %s', varargin{i});
    end
end
    
%% ====== Plot NURBS mesh ======
% --- Initializate figure ---
figure('Color','w', 'Units','normalized', 'Outerposition',[0 0 0.6 0.9]);
set(axes, 'Position', [0.1, 0.1, 0.8, 0.8]);
hold on
 
title('\bf{NURBS Surf Visualization}', 'Units','normalized', 'Position',[0.5, 1.05, 0], 'Interpreter','latex', 'fontsize', 14);
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

% --- Physical shape ---
if plot_physical_shape
    plot3(FE_grid(:,1,1), FE_grid(:,1,2), FE_grid(:,1,3), 'Color','k', 'LineWidth',3);  % edge_x_lower
    plot3(FE_grid(:,end,1), FE_grid(:,end,2), FE_grid(:,end,3), 'Color','k', 'LineWidth',3);  % edge_x_upper
    plot3(FE_grid(1,:,1), FE_grid(1,:,2), FE_grid(1,:,3), 'Color','k', 'LineWidth',3);  % edge_y_lower
    plot3(FE_grid(end,:,1), FE_grid(end,:,2), FE_grid(end,:,3), 'Color','k', 'LineWidth',3);  % edge_y_upper
    
    uistack(gca, 'top');
end

% --- Physical surf ---
if plot_physical_surf
    surf(FE_grid(:,:,1), FE_grid(:,:,2), FE_grid(:,:,3), 'FaceColor','interp', 'FaceAlpha',0.8, ...
                                                         'EdgeColor', 'interp', 'EdgeAlpha', 0.5, 'LineStyle', 'none', 'LineWidth', 1, ...
                                                         'Marker', 'none')
end

% --- Physical element ---
if plot_physical_element
    plot3(Ine(:,:,1)',Ine(:,:,2)',Ine(:,:,3)','Color','r', 'LineWidth',1.5);
end

if note_physical_element
    for i_ele = 1:size(Iee,1)
        edges = Iee(i_ele, :);
        cen_coord = [0; 0; 0];
        for i_edge = 1:length(edges)
            coord = (Ine(edges(i_edge),1,:) + Ine(edges(i_edge),end,:)) / 2;
            cen_coord = cen_coord + coord(:) / length(edges);
        end
        text(cen_coord(1), cen_coord(2), cen_coord(3), num2str(i_ele));
    end
end

if note_physical_edge
    for i_edge = 1:size(Ine,1)
        text(mean(Ine(i_edge,:,1)), mean(Ine(i_edge,:,2)), mean(Ine(i_edge,:,3)), num2str(i_edge));
    end
end



% --- Control net ---
if plot_control_net
    XYZ = reshape(CP(:,:,1:3), size(CP,1)*size(CP,1), 3);
    scatter3(XYZ(:,1), XYZ(:,2), XYZ(:,3), 'Marker', 'o', 'LineWidth',0.4, 'MarkerEdgeColor','b', 'MarkerFaceColor','b')
    
    for i_x = 1:size(CP, 1)
        plot3(CP(i_x,:,1),CP(i_x,:,2),CP(i_x,:,3), 'b--');
    end
    for i_y = 1:size(CP, 2)
        plot3(CP(:,i_y,1),CP(:,i_y,2),CP(:,i_y,3), 'b--');
    end
end

if note_control_net
    for j = 1:size(CP,2)
       for i = 1:size(CP,1) 
           text(CP(i,j,1), CP(i,j,2), CP(i,j,3), num2str((j-1)*size(CP,1)+i), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
       end
    end
end

end
