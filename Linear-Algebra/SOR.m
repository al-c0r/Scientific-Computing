clear all
%Preparing matrices as per question
A = [16, -8, -4; -8, 29, 12; -4, 12, 41];
b = [5; 25; 71];


%Determining solution
X1 = Gau_jac(A,b)
X2 = Gau_seid(A,b)
X3 = SOR(A,b,1.5) %(For SOR omega value taken as 1.5)

%Checking convergence of solution
A*X1
A*X2
A*X3

"For any w~=1, the SOR method converges slower than correspong Gauss-Seidal but faster than Gauss-Jacobi."

function sor = SOR(A,b,w)
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
            x(j) = (1-w)*x(j) + w*(b(j)-d)/A(j,j);
        end
    end
    sor = x;
end

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

function jac = Gau_jac(A,b)
    [r, r] = size(A);
    x = zeros(r,1); %Initial value
    y = 5; %No. of interations
    for i = 1:y
        c = x;
        for j = 1:r
            d = 0;
            for k = 1:r
                if k==j
                    continue
                else
                    d = d + A(j,k)*c(k);
                end
            end
            x(j) = (b(j)-d)/A(j,j);
        end
    end
    jac = x;
end