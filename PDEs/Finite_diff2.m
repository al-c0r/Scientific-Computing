clear all;

%Initialing variables
v = 50;
sig = 0.5;
k = 75;
N = 90;
M = 100;

x = zeros(1,N);
t = zeros(1,M);
a = 0; %Spatial coordinates
b = 5;
x(1) = a;
time = 10;
delta_x = (b-a)/(N-1);
delta_t = time/(N-1);

%Building x-vector 
for i=1:N-1
    x(i+1) = x(i) + (b-a)/(N-1);
end

%Building t-vector 
for i=1:M-1
    t(i+1) = t(i) + time/(N-1);
end

u = zeros(1,N);
for i = 1:N
    u(i) = exp(-x(i)^2/(2*(sig)^2))*cos(k*x(i));
end

U = zeros(N,M);
U(:,1) = u;

for j = 1:M-1 %N is spatial component and M is temporal component.
    for i = 1:N
        if i==1
            U(i,j+1) = U(i,j) - v*(delta_t)*((U(i+1,j)-U(N,j))/(2*delta_x));
        elseif i==N
            U(i,j+1) = U(i,j) - v*(delta_t)*((U(1,j)-U(i-1,j))/(2*delta_x));
        else
            U(i,j+1) = U(i,j) - v*(delta_t)*((U(i+1,j)-U(i-1,j))/(2*delta_x));
        end
    end
end

% mesh(t,x,U)
% xlabel("Distance(x) --->")
% ylabel("Time(t) --->")
% zlabel("Mass density[u(x,t)] --->")

figure(2)
v = VideoWriter("Variation of Mass density", 'MPEG-4');
v.Quality = 100;
open(v)

for i=1:M
    figure(2)
    plot(x,U(:,i))
    xlabel("Position --->")
    ylabel("Mass density --->")
    title(["Mass density at", num2str(t(i)), "Time"])
    frame = getframe(gcf);
    pause(0.2)
    writeVideo(v, frame)
end
close(v);
