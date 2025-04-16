clear all;
%Gauss elimination
A1 = [2 4 5; 1 3 5; 5 6 7]; %Solving 4c part of Assignment 1 Week2
%B1 = [23; 20; 37];

A2 = [0 -2 -3 6; -6 7 6.5 -6; 1 7.5 6.25 5.5; -12 22 15.5 -1]; %Solving 4d part of Assignment 1 Week2
%B2 = [12; -6.5; 16; 17];

%Testing the solution by multiplying it with matrix
AB1 = gau_elim(A1);
t1 = inverter(AB1)
A1*t1

AB2 = gau_elim(A2);
t2 = inverter(AB2)
A2*t2

%Converting upper triangular matrix to identity
function back = inverter(A)
    [r, c] = size(A);
    for i = 1:r
        j = i;
        while j<r
            A(r-j,:) = A(r-j,:) - A(r-j,r+1-i)*A(r+1-i,:);
            j = j+1;
        end
    end
    back = A;
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