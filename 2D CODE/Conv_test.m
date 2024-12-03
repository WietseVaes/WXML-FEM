clear all

n = 5;
Dt = .4;

NN = 6;

dx_vals = zeros(NN,1);
dt_vals = zeros(NN,1);
err = zeros(NN,1);
bnd_type = {'Kite'};

%Time_Meth = "Implicit Euler";
Time_Meth = "CN";
%Time_Meth = "BDF2";
%Time_Meth = "RK4";

for i34 = 1:NN 
    n = n*2;
    if Time_Meth == "Implicit Euler"
        Dt = Dt/4;
    elseif Time_Meth == "CN" || Time_Meth == "BDF2"
        Dt = Dt/2;
    else
        error("Method is not implemented correctly")
    end
    T = .6;

    Parameters;
    dx_vals(i34) = dx;

    [x,y,elmat,elmatbd, Id, In] = MeshShrink(bnd_type, dom_range, n, Dir_int);
    Comp
    err(i34) = L2_error(usol(:,end) -  u(:,end), elmat, x, y);
end

%Plot
figure(17);
loglog(dx_vals, err, '-o', 'LineWidth', 1.5);
xlabel('dx');
ylabel('L2 Error');
title('Error Convergence Rate');
grid on;

% Estimate convergence rate
p = polyfit(log(dx_vals), log(err), 1);
convergence_rate = p(1);

disp(['Estimated convergence rate of error: ', num2str(convergence_rate)]);


function L2_error = L2_error(err, elmat, x, y)
    
    L2_error = 0;

    
    for elem = 1:size(elmat, 1)
        vertices = elmat(elem, :);  
        x_elem = x(vertices);         
        y_elem = y(vertices);         

       
        area = abs(det([ones(3,1), x_elem, y_elem]) / 2);    
       
        u_error = (err(vertices)).^2;
        
        element_error = (area / 3) * sum(u_error); 
        L2_error = L2_error + element_error;
    end

    L2_error = sqrt(L2_error);
end

