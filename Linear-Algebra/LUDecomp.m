clear all;
%Creating matrix
A = eye(6);
for i = 1:6
    for j = 1:6
        if i==j
            A(i,j)=2;
        elseif i-j==1
            A(i,j)=-1;
        elseif j==i+1
            A(i,j)=-1;
        else
            A(i,j)=0;
        end
    end
end

B = [10 33 52 64; 23 46 75 12; 34 0 45 98; 12 43 87 56]

%U = LU(:, 1:6) %Upper triangular matrix
%L = inverter(LU(:, 7:12)) %Lower triangular matrix

U = Upper(B)
L = Lower(B)

L*U

function low = Lower(A)
    [r c] = size(A);
    B = gau_elim(A);
    low = inverter(B(:,r+1:2*r));
end

function up = Upper(A)
    [r c] = size(A);
    B = gau_elim(A);
    up = B(:,1:r);
end

%Defining inverter
function back = inverter(A)
    [r, c] = size(A);
    I = eye(r);
    AI = [A  I];
    AB = gau_elim(AI);
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
    [r, col] = size(A);
    AB = [A  eye(r)];
    [r c] = size(AB);
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
        %if AB(i,i)~=0
        %    AB(i,:) = AB(i,:)/AB(i,i);
        %end
    end
    elim = AB;
end