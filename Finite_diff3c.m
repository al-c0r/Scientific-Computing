clear all;

M = 20;
a = 0;
b = 6;

N = 30;
c = 0;
d = 8;

%Building x-vector:
x = zeros(1,M);
for i=1:M-1
    x(i+1) = x(i) + (b-a)/(M-1);
end

%Building y-vector
y = zeros(1,N);
for i=1:N-1
    y(i+1) = y(i) + (d-c)/(N-1);
end

delta_x = (b-a)/(M-1);
delta_y = (d-c)/(N-1);

f = zeros(M*N,1);
k = 1;
for i=1:M
    for j=1:N
        f(k) = x(i)^2 + y(j)^2; %sin(pi*(x(i)^2 + y(j)^2)^(1/2));
        k = k+1;
    end
end

A = zeros(M*N,M*N);
for i=1:M*N
    for j=1:M*N
        if i==j
            A(i,j) = -2*(1 + (delta_x/delta_y)^2);
        elseif mod(j-i,M*N)==1
            A(i,j) = (delta_x/delta_y)^2;
        elseif mod(j-i,M*N)==N-1
            A(i,j) = (delta_x/delta_y)^2;
        elseif mod(j-i,M*N)==N
            A(i,j) = 1;
        elseif mod(j-i,M*N)==M*N-N
            A(i,j) = 1;
        else
            A(i,j) = 0;
        end
    end
end

u = zeros(M*N,1);

A_inv = inverter(A);
u = A_inv*f;

U = zeros(M,N);

k = 1;
for i=1:M
    for j=1:N
        U(i,j) = u(k);
        k = k+1;
    end
end

mesh(x,y,U');
xlabel("Position(x)--->")
ylabel("Position(y)--->")
zlabel("Heat--->")

% Solving Linear equation:
%Converting upper triangular matrix to identity
function back = inverter(A)
    AB = gau_elim(A);
    [r, c] = size(AB);
    for i = 1:r
        j = i;
        while j<r
            AB(r-j,:) = AB(r-j,:) - AB(r-j,r+1-i)*AB(r+1-i,:);
            j = j+1;
        end
    end
    back = AB(:,r+1:2*r);
end

%Defining a Gauss-Elimination method
function elim = gau_elim(A)
    [r, column] = size(A);
    I = eye(r);
    AI = [A  I];
    [r, c] = size(AI);
    for i = 1:r
        for j = i+1:r
            if AI(i,i)==0 %Swapping if needed
                k = i+1;
                while k<r+i
                    if AI(j,i)~=0
                        d = AI(i,:);
                        AI(i,:) = AI(j,:);
                        AI(j,:) = d;
                    else
                        k = k + 1;
                    end
                end
            else
                AI(j,:) = AI(j,:) - AI(j,i)*AI(i,:)/AI(i,i); 
            end
        end
        if AI(i,i)~=0
            AI(i,:) = AI(i,:)/AI(i,i);
        end
    end
    elim = AI;
end