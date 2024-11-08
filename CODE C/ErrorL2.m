
n_values = [100,200,300]; 
dx_vals = zeros(size(n_values));
error = zeros(size(n_values));

for i=1:length(n_values)
    n = n_values(i);

    Parameters;
    dx_vals(i) = dx;

    [x,y,elmat,elmatbd, Id, In] = MeshShrink(bnd_type, dom_range, n, Dir_int);
    Comp
    error(i) = L2_error(@(x,y) usol(x,y,t(nt)), u(:,nt), elmat, x, y);
end

%Plot
figure(17);
loglog(dx_vals, error, '-o', 'LineWidth', 1.5);
xlabel('dx');
ylabel('L2 Error');
title('Error Convergence Rate');
grid on;

% Estimate convergence rate
p = polyfit(log(dx_vals), log(error), 1);
convergence_rate = p(1);

disp(['Estimated convergence rate of error: ', num2str(convergence_rate)]);


function L2_error = L2_error(u_sol_func, u_approx, elmat, x, y)
    
    L2_error = 0;

    
    for elem = 1:size(elmat, 1)
        vertices = elmat(elem, :);  
        x_elem = x(vertices);         
        y_elem = y(vertices);         

       
        area = abs(det([ones(3,1), x_elem, y_elem]) / 2);

        
        u_sol_values = u_sol_func(x_elem, y_elem);   
        
        u_approx_values = u_approx(vertices);       

        
        u_error = (u_sol_values - u_approx_values).^2;
        
        element_error = (area / 3) * sum(u_error); 
        L2_error = L2_error + element_error;
    end

    L2_error = sqrt(L2_error);
end

