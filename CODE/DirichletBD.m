function DB_tf = DirichletBD(bd_x,bd_y,s)
% Dirichlet Boundary condition is defined here.
% Gives back 1 if (x,y) is on the Dirichlet boundary condition
%            0 if not.
% x and y can be a vector, if so, tf is a vector of 1 and 0s as described
%above
DB_int = linspace(0,0.25*pi,length(bd_x));
DB_points = s.Zorigional(DB_int);
DB_tf = zeros(size(bd_x));

tol = 1e-6;

for i = 1:length(bd_x)
    difference = abs(DB_points - (bd_x(i) + 1i*bd_y(i)));
    if any(difference < tol)
        DB_tf(i) = 1;
    end
end

end