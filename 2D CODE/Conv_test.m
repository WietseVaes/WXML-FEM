clear all

n = 5;
Dt = .4;

NN = 5;

dx_vals = zeros(NN,1);
dt_vals = zeros(NN,1);
err = zeros(NN,1);
bnd_type = {'Kite'};

% Time_Meth = "Implicit Euler";
% Time_Meth = "CN";
%Time_Meth = "BDF2";
%Time_Meth = "RK4";
Time_Meth = ["Implicit Euler", "CN", "BDF2"];
counter = 1;
all_err = zeros(3,NN);
all_dx_vals = zeros(3,NN);

for met = Time_Meth
    dx_vals = zeros(NN,1);
    dt_vals = zeros(NN,1);
    err = zeros(NN,1);
    n = 5;
    Dt = .4;
    for i34 = 1:NN 
        n = n*2;
        if met == "Implicit Euler"
            Dt = Dt/4;
        elseif met == "CN" || met == "BDF2"
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
    all_err(counter,:) = err;
    all_dx_vals(counter,:) = dx_vals;
    counter = counter +1;
end
save('error_and_dx_results.mat', 'all_err','all_dx_vals');

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

