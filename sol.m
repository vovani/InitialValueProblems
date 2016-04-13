clear;
t = [0, 1];
x = [0, 1];
A = 1i .* [0 1; 1 0];
C = 1i .* [3, -1; -1, 3];
w = 2 * pi;

boundary_cond = @(x) [cos(w .* x); - 1i .* sin(w .* x)];
pde_sol = @(x, t) expm((C - (w ^ 2) .* A) * t) * boundary_cond(x);
dt = @(pow, lambda) @(dx) (dx ^ pow) * lambda;

% Ns = 2 .^ (4:7);
% plot_errors('FE : dt = 0.1 * (dx^3)', Ns, dt(3, 0.5), x, t, @(some_grid) fe(some_grid, A, C, boundary_cond), pde_sol);

% Ns = 2 .^ (4:7);
% plot_errors('FE : dt = 1 * (dx^4)', Ns, dt(4, 1), x, t, @(some_grid) fe(some_grid, A, C, boundary_cond), pde_sol);

% figure; 
% Ns = 2 .^ (7:12);
% plot_errors('BE : dt = lam * dx', Ns, dt(1,  0.1), x, t, @(some_grid) be(some_grid, A, C, boundary_cond), pde_sol);

% figure;
% Ns = 2 .^ (4:8);
% plot_errors('BE : dt = 0.1 * dx^2', Ns, dt(2, 0.1), x, t, @(some_grid) df(some_grid, A, C, boundary_cond), pde_sol);
% 
% figure;
% Ns = 2 .^ (4:8);
% plot_errors('CN : dt = 0.1 * dx', Ns, dt(1, 0.1), x, t, @(some_grid) cn(some_grid, A, C, boundary_cond), pde_sol);
% 
% figure;
% Ns = 2 .^ (4:8);
% plot_errors('CN : dt = dx^2', Ns, dt(2, 1), x, t, @(some_grid) cn(some_grid, A, C, boundary_cond), pde_sol);
% 
% figure;
% Ns = 2 .^ (4:8);
% plot_errors('LF : dt = 0.2 * dx ^ 2', Ns, dt(2, 0.2), x, t, @(some_grid) lf(some_grid, A, C, boundary_cond), pde_sol);
% 
% figure;
% Ns = 2 .^ (4:8);
% plot_errors('LF : dt = 0.25 * dx ^ 2', Ns, dt(2, 0.25), x, t, @(some_grid) lf(some_grid, A, C, boundary_cond), pde_sol);
% % 
figure;
Ns = 2 .^ (4:8);
plot_errors('DF : dt = 0.1 * dx ^ 2', Ns, dt(1 , 0.1), x, t, @(some_grid) df(some_grid, A, C, boundary_cond), pde_sol);
% figure;
% Ns = 2 .^ (4:8);
% plot_errors('DF : dt = 0.1 * dx ^ 2', Ns, dt(2 , 0.1), x, t, @(some_grid) df(some_grid, A, C, boundary_cond), pde_sol);