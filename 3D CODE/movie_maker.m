function [] = movie_maker(u, x, y, z, sp, file_name)

nt = size(u,2);
n = round((length(x))^(1/3));

% For irregular shape we will need to worry about a keep filter, but we'll
% get there when we get there.

for i1 = sp:nt
    %for room make all x,y,z, the original versions
    uplt = u(:,i1);

    xx = reshape(x,n,n,n); yy = reshape(y,n,n,n); zz = reshape(z,n,n,n);

    pltvals = reshape(uplt,n,n,n);

    vtk_name =['3D CODE/VTK/',file_name,'_',num2str(i1),'.vtk'];
    vtkwrite(vtk_name,'structured_grid',xx,yy,zz,'scalars','Wave',pltvals);
end

end

% function [] = movie_maker(u, x, y, z, sp, file_name)
% 
% nt = size(u, 2);
% n = round((length(x))^(1/3));
% 
% 
% keep = true(size(x)); % Initialize
% 
% for i = 1:length(x)
%     if y(i) > 2 * x(i) + 5
%         keep(i) = false;
%     end
% 
%     if y(i) > -2 * x(i) + 13
%         keep(i) = false;
%     end
% end
% 
% x = x(keep);
% y = y(keep);
% z = z(keep);
% u = u(keep, :);
% 
% n_new = round((length(x))^(1/3));
% size(x)
% for i1 = sp:nt
%     vals = u(:, i1);
% 
%     % Reshape using the new dimensions
%     xx = reshape(x, n_new, n_new, n_new);
%     yy = reshape(y, n_new, n_new, n_new);
%     zz = reshape(z, n_new, n_new, n_new);
%     pltvals = reshape(vals, n_new, n_new, n_new);
% 
%     % Save the file
%     vtk_name = ['3D CODE/VTK/', file_name, '_', num2str(i1), '.vtk'];
%     vtkwrite(vtk_name, 'structured_grid', xx, yy, zz, 'scalars', 'Wave', pltvals);
% end
% 
% end
% 
% 
