function [v] = be(a_grid, A, C, boundery_func)
    N = a_grid.N;
 
    DpDn = sparse(gallery('circul', [-2 1 zeros(1, N - 3) 1])) ./ (a_grid.dx ^ 2);

    Q = speye(2 * N) - a_grid.dt .* (kron(A, DpDn) + kron(C, speye(N)));
    % Q = sparse(eye(2 * N) - a_grid.dt .* kron(A, DpDn)));
    
    v = boundery_func(a_grid.x);
    v = [v(1, :), v(2, :)].';
    %v = (inv(Q) ^ size(a_grid.t, 2)) * v;
    for j = 1:size(a_grid.t, 2)
       v = Q \ v; 
    end
    % v = ((inv(Q) * (eye(2*N) + a_grid.dt .* kron(C, eye(N)))) ^ size(a_grid.t, 2)) * v;
    v = [v(1:N, 1), v((N + 1):2 * N, 1)].';
end