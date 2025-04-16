clear all;

%Creating matrix
A = [4 -2 1; 3 6 -1; 3 5 9];
b = [5; 25; 71];

U = Upper(A)
L = Lower(A)
D = Diag(A)

L*U %Check if their product is actually the matrix itself or not

%Defining Lower & Upper matrix converter
function low = Lower(A)
    [r c] = size(A);
    B = gau_elim(A);
    low = inverter(B(:,r+1:2*r));
end

function diag = Diag(A)
    [r, c] = size(A);
    B = eye(r);
    C = gau_elim(A);
    for i = 1:r
        B(i,i)=C(i,i);
    end
    diag = B
end

function up = Upper(A)
    [r c] = size(A);
    B = gau_elim(A);
    up = B(:,1:r);
    for i=1:r
        B(:,i) = B(:,i)/B(i,i)
end

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
        %if AI(i,i)~=0
        %   AI(i,:) = AI(i,:)/AI(i,i);
        %end
    end
    elim = AI;
end