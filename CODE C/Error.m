%Calculate Error 

n_values = [200, 400, 600]; 
h_vals = zeros(size(n_values));
max_error = zeros(size(n_values)); 
avg_error = zeros(size(n_values));

% Loop over each n
for i = 1:length(n_values)
    
    
    Parameters;
    n = n_values(i);
  
    [x,y,elmat,elmatbd, Id, In] = MeshShrink(bnd_type, dom_range, n, Dir_int,f,g,h); %edit the inputs

    
    Comp
    
    dx = abs((dom_range{1}(1) - dom_range{1}(2))/n);
    dy = abs((dom_range{2}(1) - dom_range{2}(2))/n);
    
    h = max(dx,dy);
    h(i) = h;
    e = abs(usol(x,y, t(nt))-u(:,t(nt)));

    avg_e = meanabs((e)) * h; %l1 norm
    max_e = max(e);

    max_error(i) = max_e;
    avg_error(i) = avg_e;
end

%max
figure(4);
loglog(h_vals, max_error, '-o', 'LineWidth', 1.5);
xlabel('Mesh size');
ylabel('Error');
title('Error Convergence Rate');
grid on;

% Estimate convergence rate
p = polyfit(log(h_vals), log(max_error), 1);
max_convergence_rate = p(1);

disp(['Estimated convergence rate of max error: ', num2str(max_convergence_rate)]);


%avg
figure(5);
loglog(h_vals, avg_error, '-o', 'LineWidth', 1.5);
xlabel('Mesh size');
ylabel('Error');
title('Error Convergence Rate');
grid on;


p = polyfit(log(h_vals), log(avg_error), 1);
avg_convergence_rate = p(1);

disp(['Estimated convergence rate of avg error: ', num2str(avg_convergence_rate)]);
