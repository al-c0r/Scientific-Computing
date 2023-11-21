clear all;

N = 20;
a = 0;
b = 6;

c = 0;
d = 6;

%Building x-vector:
x = zeros(1,N);
for i=1:N-1
    x(i+1) = x(i) + (b-a)/(N-1);
end

%Building y-vector
y = zeros(1,N);
for i=1:N-1
    y(i+1) = y(i) + (d-c)/(N-1);
end

delta_x = (b-a)/(N-1);

f = zeros(N^2,1);
k = 1;
for i=1:N
    for j=1:N
        f(k) = sin(pi*(x(i)^2 + y(j)^2)^(1/2));
        k = k+1;
    end
end

A = zeros(N^2,N^2);
for i=1:N^2
    for j=1:N^2
        if i==j
            A(i,j) = -4;
        elseif mod(j-i,N^2)==1
            A(i,j) = 1;
        elseif mod(j-i,N^2)==N-1
            A(i,j) = 1;
        elseif mod(j-i,N^2)==N
            A(i,j) = 1;
        elseif mod(j-i,N^2)==N^2-N
            A(i,j) = 1;
        else
            A(i,j) = 0;
        end
    end
end

u = zeros(N^2,1);
size(x)
size(y)
size(u)


A_inv = inverter(A);
u = A_inv*b;
Error = A*u - b;
U = zeros(N,N);
k = 1;
for i=1:N
    for j=1:N
        U(i,j) = u(k);
        k = k+1;
    end
    k = k+1;
end

mesh(x,y,U);
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