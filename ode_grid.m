classdef ode_grid
   properties
      N
      x
      dx
      dt
      t
   end
   methods
      function obj = ode_grid(x, t, N, dt_calc)
         obj.N = N;
         obj.x = linspace(x(1), x(2), N);
         obj.dx = (x(2) - x(1)) / (N - 1);
         obj.dt = dt_calc(obj.dx);
         obj.t = (t(1) + obj.dt) : obj.dt : t(2);
         if obj.tf() ~= t(2)
             error('time interval is not accurate!');
         end
      end
      function r = tf(obj)
         r = obj.t(length(obj.t)); 
      end
      function r = t0(obj)
         r = obj.t(0); 
      end
   end
end