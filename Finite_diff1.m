clear all;
%Initialising variables
alp = 0.62; %alpha = k(delta(t))/(delta(x)^2)
N = 150; %No. of space parts
M = 100; %No. of time parts

x = zeros(1,N);
t = zeros(1,M);
a = 0; %Spatial coordinates
b = 5;
x(1) = a;
time = 10;

%Building x-vector 
for i=1:N-1
    x(i+1) = x(i) + (b-a)/(N-1);
end

%Building t-vector 
for i=1:M-1
    t(i+1) = t(i) + time/(N-1);
end

u = sin(pi*x);
U = zeros(M,N); %M-rows (temporal-component) and N-columns (spatial component)
U(M,:) = u;
for j=1:M-1
    for i=2:N-1
        U(M-j,i) = U(M-j+1,i) + alp*(U(M-j+1,i+1) + U(M-j+1,i-1) - 2*U(M-j+1,i));
    end
end

mesh(x,t,U)
xlabel("Distance(x) --->")
ylabel("Time(t) --->")
zlabel("Heat[u(x,t)] --->")