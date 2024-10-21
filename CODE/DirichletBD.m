function DB_tf = DirichletBD(bd_x, bd_y, s)
    % DirichletBD: Determines if boundary points are on the Dirichlet boundary
    % based on the parameterization 't' of the boundary curve.
    %
    % bd_x, bd_y: Boundary point x and y coordinates
    % s: Structure containing the original boundary parameterization
    % Returns:
    % DB_tf: 1 for Dirichlet boundary points (0 ≤ t ≤ π), 0 for Neumann points (π < t ≤ 2π)

    DB_tf = zeros(size(bd_x));
    tol = 1e-6;
    
    for i = 1:length(bd_x)
        % Find the corresponding 't' value for the boundary point (bd_x(i), bd_y(i))
        point = bd_x(i) + 1i * bd_y(i);  % Complex representation of the boundary point
        difference = abs(s.original(s.t) - point);  % Difference with boundary points
        
        % Find the closest parameter value in 't' that matches the point
        [~, idx] = min(difference);
        t_value = s.t(idx);  % The 't' value closest to the boundary point
        
        % Check if the 't' value is in the Dirichlet range [0, pi]
        if t_value >= 0 && t_value <= pi
            DB_tf(i) = 1;  % This point is on the Dirichlet boundary
        end
    end
end