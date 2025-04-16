clear all

A = [1, -1; 1, 1; -1, -1]
Q = Gram_Sch(A)

R = inv_up(A)

Q*R

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