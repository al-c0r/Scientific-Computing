clear all;
%Initialising variables
k = 0.62; %alpha = k(delta(t))/(delta(x)^2)
N1 = 50; %No. of space parts (X)
N2 = 60; %(Y)
N3 = 70; %(Z)
M = 80; %No. of time parts

x = zeros(1,N1);
y = zeros(1,N2);
z = zeros(1,N3);
t = zeros(1,M);

a1 = 0; %Spatial coordinates (X)
a2 = 5; 
b1 = 0; %(Y)
b2 = 6;
c1 = 0; %(Z)
c2 = 7;
time = 10;

delta_x = (a2-a1)/(N1-1);
delta_y = (b2-b1)/(N2-1);
delta_z = (c2-c1)/(N3-1);
delta_t = time/(M-1);


x(1) = a1;
y(1) = b1;
z(1) = c1;

%Building x-vector 
for i=1:N1-1
    x(i+1) = x(i) + (a2-a1)/(N1-1);
end

%Building y-vector 
for i=1:N2-1
    y(i+1) = y(i) + (b2-b1)/(N2-1);
end

%Building z-vector 
for i=1:N3-1
    z(i+1) = z(i) + (c2-c1)/(N3-1);
end


%Building t-vector 
for i=1:M-1
    t(i+1) = t(i) + time/(M-1);
end

u = zeros(N1,N2,N3); %+ y^2 + z^2)^(1/2));

for i=1:N1
    for j=1:N2
        for k=1:N3
            u(i,j,k) = sin(pi*(x(i)^2 + y(j)^2 + z(k)^2)^(1/2));
        end
    end
end

U = zeros(N1,N2,N3,M); %M-rows (temporal-component) and N-columns (spatial component)
U(:,:,:,M) = u;

for l=1:M-1
    for i=2:N1-1
        for j=2:N2-1
            for k=2:N3-1
                U(i,j,k,M-l) = (1 - 2*k*(delta_t)*(1/(delta_x)^2 + 1/(delta_y)^2 + 1/(delta_z)^2))*U(i,j,k,M-l+1) + k*(delta_t)*(((U(i+1,j,k,M-l+1)+U(i-1,j,k,M-l+1))/(delta_x)^2) + (U(i,j+1,k,M-l+1)+U(i,j-1,k,M-l+1))/(delta_y)^2 + (U(i,j,k+1,M-l+1)+U(i,j,k-1,M-l+1))/(delta_z)^2);
            end
        end
    end
end

Ux = zeros(N1,M); %For plotting U with X and T while fixing Y and Z
for i=1:N1
    for j=1:M
        Ux(i,j)=U(i,N2,N3,j);
    end
end

Uy = zeros(N2,M); %For plotting U with Y and T while fixing Z and X
for i=1:N2
    for j=1:M
        Uy(i,j)=U(N1,i,N3,j);
    end
end

Uz = zeros(N3,M); %$For plotting U with Z and T while fixing X and Y
for i=1:N3
    for j=1:M
        Uz(i,j)=U(N1,N2,i,j);
    end
end

figure(1)
mesh(t,x,Ux)
ylabel("Distance(x) --->")
xlabel("Time(t) --->")
zlabel("Heat[U(x,t)] (y=N2,z=N3) --->")
hold on

figure(2)
mesh(t,y,Uy)
ylabel("Distance(y) --->")
xlabel("Time(t) --->")
zlabel("Heat[u(y,t)] (z=N3,x=N1)--->")
hold on
 
figure(3)
mesh(t,z,Uz)
ylabel("Distance(z) --->")
xlabel("Time(t) --->")
zlabel("Heat[U(z,t)] (x=N1,y=N2) --->")
