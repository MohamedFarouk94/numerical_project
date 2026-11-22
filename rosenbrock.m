function f = rosenbrock(x)
%ROSENBROCK  Rosenbrock "banana" function.
%
%   f = ROSENBROCK(x) evaluates
%
%       f(x1, x2) = (1 - x1)^2 + 100 * (x2 - x1^2)^2
%
%   at the point x = [x1; x2]. The global minimum is f(1, 1) = 0.

    f = (1 - x(1))^2 + 100 * (x(2) - x(1)^2)^2;
end
