clear all;

A = [0.5, 0, 0.2; 0, 3.15, -1; 0.57, 0, -7.43];
QR_eig(A)
"Therefore, the eigenvalues by QR-Algorithm is"
[A(1,1), A(2,2), A(3,3)]

function eig = QR_eig(A)
    n = 100; %No. of iterations
    [m,m] = size(A);
    B = zeros(m,m);
    B = Gram_Sch(A)*inv_up(A);
    for i = 1:10
        B = inv_up(B)*Gram_Sch(B);
    end
    eig = B;
end     

function R = inv_up(A)
    [n, n] = size(A);
    C = zeros(n,n);
    G = Gram_Sch(A);
    for i = 1:n
        for j = 1:n
            if i>j
                C(i,j) = 0;
            else
                C(i,j) = A(:,j)'*G(:,i);
            end
        end
    end
    R = C;
end    

function gram = Gram_Sch(A)
    [m, n] = size(A);
    B = zeros(m,n);
    for i = 1:n
        x = zeros(m,1);
        for j = 1:i
            %y = A(:,j);
            if j == i
                x = x +  A(:,j);
            else
                x = x - (A(:,i)'*B(:,j))*B(:,j);
            end
        end
        B(:,i) = x/(sqrt(x'*x));
    end
    gram = B;
end