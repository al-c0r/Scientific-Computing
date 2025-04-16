clear all
%Preparing matrices as per question
A1 = [4, -2, 1; 3, 6, -1; 3, 5, 9];
b1 = [5; 25; 71];

A2 = [2, 1, 3; 3, 4 -2; 4, 1, 1];
b2 = [7; 9; 7];

A3 = [10, -1, 2, 0; -1, 11, -1, 3; 2, -1, 10, -1; 0, 3, -1, 8;];
b3 = [6; 25; -11; 15];

%Determining solution
X1 = Gau_seid(A1,b1)
X2 = Gau_seid(A2,b2)
X3 = Gau_seid(A3,b3)

%Checking convergence of solution
A1*X1
A2*X2 %This one doesn't converge at all!
A3*X3

function jac = Gau_seid(A,b)
    [r, r] = size(A);
    x = zeros(r,1); %Initial value
    y = 5; %No. of interations
    for i = 1:y
        for j = 1:r
            d = 0;
            for k = 1:r
                if k==j
                    continue
                else
                    d = d + A(j,k)*x(k);
                end
            end
            x(j) = (b(j)-d)/A(j,j);
        end
    end
    jac = x;
end

