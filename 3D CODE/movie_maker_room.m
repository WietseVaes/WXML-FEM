function [] = movie_maker_room(u, x, y, z, sp, file_name, keep)
    % Ensure irregular shapes work by handling unstructured grids for VTK.
    
    nt = size(u, 2); % Number of timesteps

    % Iterate over timesteps
    for i1 = sp:nt
        vals = u(:, i1);

        % Create the VTK file for this timestep
        vtk_name = ['3D CODE/VTK/', file_name, '_', num2str(i1), '.vtk'];
        
        % Use unstructured grid for irregular shapes
        %write_irregular_vtk(vtk_name, x, y, z, vals);
        vtkwrite(vtk_name,'unstructured_grid',x,y,z,'scalars','Wave',vals,vals,vals);

    end
end

% function write_irregular_vtk(file_name, x, y, z, connectivity, scalars)
%     % Write a VTK file for an unstructured grid with non-symmetric shapes.
%     % 
%     % points: Nx3 array of [x, y, z] coordinates
%     % connectivity: MxK array of indices defining cells (triangles, tetrahedra, etc.)
%     % scalars: Nx1 array of scalar values (e.g., wave amplitudes) at each point
% 
%     n_points = size(x);
%     n_cells = size(connectivity, 1);
% 
%     % Open file
%     fid = fopen(file_name, 'w');
% 
%     % Write VTK Header
%     fprintf(fid, '# vtk DataFile Version 2.0\n');
%     fprintf(fid, 'Unstructured Grid Example\n');
%     fprintf(fid, 'ASCII\n\n');
%     fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
% 
%     % Write Points
%     fprintf(fid, 'POINTS %d float\n', n_points);
%     for i = 1:n_points
%         fprintf(fid, '%f %f %f\n', x(i), y(i), z(i));
%     end
% 
%     % Write Connectivity (Cells)
%     fprintf(fid, '\nCELLS %d %d\n', n_cells, n_cells * (size(connectivity, 2) + 1));
%     for i = 1:n_cells
%         fprintf(fid, '%d ', size(connectivity, 2)); % Number of points in the cell
%         fprintf(fid, '%d ', connectivity(i, :) - 1); % Indices (0-based for VTK)
%         fprintf(fid, '\n');
%     end
% 
%     % Write Cell Types
%     fprintf(fid, '\nCELL_TYPES %d\n', n_cells);
%     cell_type = 5; % Example: triangles (VTK_TRIANGLE = 5)
%     fprintf(fid, '%d\n', cell_type * ones(1, n_cells));
% 
%     % Write Scalars
%     fprintf(fid, '\nPOINT_DATA %d\n', n_points);
%     fprintf(fid, 'SCALARS Wave float 1\n');
%     fprintf(fid, 'LOOKUP_TABLE default\n');
%     for i = 1:n_points
%         fprintf(fid, '%f\n', scalars(i));
%     end
% 
%     % Close file
%     fclose(fid);
% end
% 
