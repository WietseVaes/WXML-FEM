function [] = movie_maker(u, x, y, z, sp, file_name)

nt = size(u,2);
n = round((length(x))^(1/3));

% For irregular shape we will need to worry about a keep filter, but we'll
% get there when we get there.

for i1 = sp:nt

    vals = u(:,i1);
    
    xx = reshape(x,n,n,n); yy = reshape(y,n,n,n); zz = reshape(z,n,n,n);

    pltvals = reshape(vals,n,n,n);
    
    vtk_name =['3D CODE/VTK/',file_name,'_',num2str(i1),'.vtk'];
    vtkwrite(vtk_name,'structured_grid',xx,yy,zz,'scalars','Wave',pltvals);
end

end