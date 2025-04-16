clear all;

A = [4, 1, 1; 1, 3, 1; 1, 1, 5];

B = [1, 0, 0; 0, 1, 0; 0, 0, 1]; %3 vectors that has to A-orthonormalise

"The A-conjugate orthonormalised vectors are"

C = A_Gram_Sch(B,A)

function gram = A_Gram_Sch(B,A)
    [m, n] = size(B);
    D = zeros(m,n);
    for i = 1:n
        x = zeros(m,1);
        for j = 1:i
            if j == i
                x = x +  B(:,j);
            else
                x = x - Inner_A_product(B(:,i),A,D(:,j))*D(:,j);
            end
        end
        D(:,i) = x/(sqrt(Inner_A_product(x,A,x)));
    end
    gram = D;
end

function in_Ap = Inner_A_product(x,A,y) %Defining A-conjugate product
    in_Ap = x'*A*y;
end
