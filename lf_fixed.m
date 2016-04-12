function [v] = lf_fixed(a_grid, A, C, boundery_func)
    N = a_grid.N;
     
    DpDn = sparse(gallery('circul', [-2 1 zeros(1, N - 3) 1])) ./ (a_grid.dx ^ 2);
    Ep = sparse(gallery('circul', [0 1 zeros(1, N - 2)]));
    En = sparse(gallery('circul', [zeros(1, N - 1) 1]));
    FE_Q = sparse(eye(2*N) + a_grid.dt .* (kron(A, DpDn) + kron(C, eye(N))));

    Q11 = 2 * a_grid.dt .* kron(A, DpDn) + a_grid.dt .* kron(C, Ep + En);
    Q = [Q11 eye(2 * N); eye(2 * N) zeros(2 * N)];
    
    v1 = boundery_func(a_grid.x);
    v1 = [v1(1, :), v1(2, :)].';
    v2 = FE_Q * v1;
    v = [v2(:, 1); v1(:, 1)];
    
    v = (Q ^ (size(a_grid.t, 2) - 1)) * v;
    v = horzcat(v(1:N, 1), v((N + 1):2 * N, 1)).';
end