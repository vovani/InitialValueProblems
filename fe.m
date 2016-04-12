function [v] = fe(a_grid, A, C, boundery_func)
    N = a_grid.N;
      
    DpDn = sparse(gallery('circul', [-2 1 zeros(1, N - 3) 1])) ./ (a_grid.dx ^ 2);
    Q = sparse(eye(2*N) + a_grid.dt .* (kron(A, DpDn) + kron(C, eye(N))));
    v = boundery_func(a_grid.x);
    v = [v(1, :), v(2, :)].';
    v = (Q ^ size(a_grid.t, 2)) * v;
    v = horzcat(v(1:N, 1), v((N + 1):2 * N, 1)).';
end