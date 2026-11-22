function H = hessian(x)
%HESSIAN  Analytical Hessian of the Rosenbrock function.
%
%   H = HESSIAN(x) returns the symmetric 2x2 Hessian matrix at
%   x = [x1; x2]:
%
%       [ 2 - 400 x2 + 1200 x1^2,   -400 x1 ]
%       [ -400 x1,                   200     ]

    x1 = x(1);
    x2 = x(2);

    H = [2 - 400 * x2 + 1200 * x1^2, -400 * x1;
         -400 * x1,                   200];
end
