function g = grad(x)
%GRAD  Analytical gradient of the Rosenbrock function.
%
%   g = GRAD(x) returns the 2x1 gradient vector at x = [x1; x2]:
%
%       df/dx1 = -2 (1 - x1) - 400 x1 (x2 - x1^2)
%       df/dx2 =  200 (x2 - x1^2)

    x1 = x(1);
    x2 = x(2);

    g = [-2 * (1 - x1) - 400 * x1 * (x2 - x1^2);
          200 * (x2 - x1^2)];
end
