function [error] = calc_approx_error(u, v, dx)
% u - expected result vector, u(*, 1) is the result at leftmost point of
%     the grid.
% v - actual result vector, same format as u
% dx - commonly reffered as "h".
    error = sqrt(dx * sum( sum (abs(u - v).^2) ));
end