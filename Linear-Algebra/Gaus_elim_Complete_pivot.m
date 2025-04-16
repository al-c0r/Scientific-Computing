%Gauss elimination
A1 = [2 4 5; 1 3 5; 5 6 7]; %Solving 4c part of Assignment 1 Week2
B1 = [23; 20; 37];

A2 = [0 -2 -3 6; -6 7 6.5 -6; 1 7.5 6.25 5.5; -12 22 15.5 -1]; %Solving 4d part of Assignment 1 Week2
B2 = [12; -6.5; 16; 17];

%Testing the solution by multiplying it with matrix
AB1 = gau_elim(A1,B1)
t1 = back_sub(AB1)
A1*t1'

AB2 = gau_elim(A2,B2)
t2 = back_sub(AB2)
A2*t2'

%Defining Backsubstition
function back = back_sub(AB)
    [r c] = size(AB)
    for i = 1:r
        y = AB(r+1-i , c);
        j = 1;
        while j<i
            y = y - AB(r+1-i, c-j)*x(r+1-j);
            j = j+1;
        end
        x(r+1-i) = y;
    end
    back = x;
end

%Defining a Gauss-Elimination method
function elim = gau_elim(A,B)
    AB = [A  B];
    [r, c] = size(AB);
    for i = 1:r %Partial Pivoting if required.
        if abs(AB(1,1))<abs(AB(i,1))
            f = AB(1,:);
            AB(1,:) = AB(i,:);
            AB(i,:) = f;
        end
    end

    for i = 1:r
        AB(i,:) = AB(i,:)/maximum(AB(i,:)); %Normalising by dividing it with it's max.
    end

    for i = 1:r
        for j = i+1:r
            if AB(i,i)==0 %Swapping if needed
                k = i+1;
                while k<r+i
                    if AB(j,i)~=0
                        d = AB(i,:);
                        AB(i,:) = AB(j,:);
                        AB(j,:) = d;
                    else
                        k = k + 1;
                    end
                end
            else
                AB(j,:) = AB(j,:) - AB(j,i)*AB(i,:)/AB(i,i); 
            end
        end
        if AB(i,i)~=0
            AB(i,:) = AB(i,:)/AB(i,i);
        end
    end
    elim = AB;
end

function max = maximum(x)
    g = x(1);
    for i = 1:length(x)
        if abs(x(i)) > abs(g)
            g = x(i);
        end
    end
    max = g;
end