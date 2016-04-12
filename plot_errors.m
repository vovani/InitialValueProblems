function plot_errors(name, Ns, calc_dt, x, t, method, pde_sol)
    fprintf('\n\n')
    disp(['Showing solution for ', name])
    
    errors = zeros(size(Ns));
    for j = 1 : length(Ns)
        N = Ns(j);
        a_grid = ode_grid([x(1), x(2) - 1/N], t, N, calc_dt);
        v = method(a_grid);
        u = pde_sol(a_grid.x, a_grid.tf());
        errors(j) = calc_approx_error(u, v, a_grid.dx);
        fprintf('N = %d, dt = %f, tf = %d, err = %d\n', ... 
                    N, a_grid.dt, a_grid.tf, errors(j));
    end
    
    a_fit= polyfit(log2(1 ./ Ns), log2(errors), 1); 
    alpha = a_fit(1);
    beta = a_fit(2);
    hold on;
    xx = linspace(log2(1 / Ns(1)), log2(1 / Ns(size(Ns, 2))), 1000);
    plot(log2(1 ./ Ns), log2(errors), 'b*', xx, beta + alpha .* xx, 'r-');
    hold off;
    str = strcat('alpha = ', num2str(alpha));
    title([name , ', ' , str]);
    xlabel('log(dx)');
    ylabel('log(||E||)');
end