clear all;

A = [0.5, 0, 0.2; 0, 3.15, -1; 0.57, 0, -7.43];
Dom_eigval(A)

function pow = Dom_eigval(A)
    n = 100; %No. of iterations
    x = [1;1;1];
    for i=1:n
        x = A*x;
    end
    pow = x'*A*x/(x'*x);
end
