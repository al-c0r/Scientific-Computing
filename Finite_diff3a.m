clear all;

N = 20;
a = 0;
b = 6;

%Building x-vector:
x = zeros(1,N);
for i=1:N-1
    x(i+1) = x(i) + (b-a)/(N-1);
end

delta_x = (b-a)/(N-1);

f = sin(pi*x);

A = zeros(N,N);
for i=1:N
    for j=1:N
        if i==j
            A(i,j) = -2;
        elseif abs(i-j)==1
            A(i,j) = 1;
        else
            A(i,j) = 0;
        end
    end
end

size(A)
b = zeros(N,1);
u = zeros(N,1);

%Bounds
f1 = 0;
f2 = 10;
for i=1:N
    if i==1
        b(i) = ((delta_x)^2)*f(i) - f1;
    elseif i==N
        b(i) = ((delta_x)^2)*f(i) - f2;
    else
        b(i) = ((delta_x)^2)*f(i);
    end
end

A_inv = inverter(A);
u = A_inv*b;
Error = A*u - b

plot(x,u);
xlabel("Position(x)--->")
ylabel("Heat--->")

%Solving Linear equation:
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