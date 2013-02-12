% -- Modifiable Parameters -- %
m = 5;                      % Mass
c = [0,0,0;0,0,0;0,0,0];    % First Moment of Inertia Matrix
J = [10,0,0;0,10,0;0,0,20]; % Second Moment of Inertia
v0_lin = [1 -2 0]';          % Initial Linear Velocity
v0_ang = [0 1 -5]';     % Initial Angular Velocity
v0 = [ v0_lin; v0_ang ];    % Augmented Initial Velocity Vector
f_ext = @(t) [0 t 0]';      % External Force Function
g_ext = @(t) [-2*t 0 1]';      % External Torque Function

% -- Computed Parameters -- %
M = [ m*eye(3) -1*skew_sym(c); skew_sym(c) J ];     % System Inertia Matrix
p0 = M*v0;                                          % Compute Intial Momentum

% -- Numerical Simultaion Parameters -- %
a = 0;                                      % Start integration region
b = 10;                                     % End integration region
multiplier = 10;                           % Factor to multiply ideal step count
stepFunc = @(a,b,v0) 2*norm(v0)*(b-a);      % Function to compute number of steps
n = round(multiplier*stepFunc(a,b,v0_ang)); % Compute Number of Steps
tolerance = 0.00000001;                            % Tolerance for adaptive schemes
%n = 10 ;                                  % Explicit number of steps

% -- Integration Procedure -- %
% Types:    1 - RK4     | Classical Runge-Kutta 4th Order Method
%           2 - IE      | Implicit (Backwards) 1st Order Method
%           3 - GRK4    | Implicit Gauss-Legendre 4th Order Method
%           4 - DOP54   | Adaptive Dormand-Prince 5/4th Order Method
type = 1;                                   % Method Type
if type == 1     %RK4
    [R,P,V,T] = rk4_freespace(a,b,p0,v0,f_ext,g_ext,M,n);
elseif type == 2 %IE
    [R,P,V,T] = impliciteuler_freespace(a,b,p0,v0,f_ext,g_ext,M,n);
elseif type == 3 %GRK4
    [R,P,V,T] = grk4_freespace(a,b,p0,v0,f_ext,g_ext,M,n);
elseif type == 4 %DOP54
    [R,P,V,T,n] = dopri54_freespace(a,b,p0,v0,f_ext,g_ext,M,tolerance);
elseif type == 5 %AGRK4
    [R,P,V,T,n] = agrk4_freespace(a,b,p0,v0,f_ext,g_ext,M,tolerance);
elseif type == 6 %RK
    [R,P,V,T] = rk_freespace(a,b,p0,v0,f_ext,g_ext,M,n,'Gauss 2');
end
% -- Format and Display Results -- %
fileID = fopen('results_RK.txt','wt');
fprintf(fileID,' Time   Momentum     Norm(P)  Velocity    Norm(V)       Type\n\n');
for i=1: 1:n+1
    fprintf(fileID,'%6.3f%10.3g\t    %10.3g \t\t\tLinear\n',R(i,1),R(i,2),R(i,8));
    fprintf(fileID,'      %10.3g\t    %10.3g \t\t\tLinear\n',R(i,3),R(i,9));
    fprintf(fileID,'      %10.3g %10.3g %10.3g %10.3g \tLinear\n\n',R(i,4),norm(R(i,2:4)),R(i,10),norm(R(i,8:10)));
    fprintf(fileID,'      %10.3g\t    %10.3g \t\t\tAngular\n',R(i,5),R(i,11));
    fprintf(fileID,'      %10.3g\t    %10.3g \t\t\tAngular\n',R(i,6),R(i,12));
    fprintf(fileID,'      %10.3g %10.3g %10.3g %10.3g \tAngular\n\n',R(i,7),norm(R(i,5:6)),R(i,13),norm(R(i,11:13)));
end
Vlin_norm = zeros(1,n+1);
for i=1:n+1
    Vlin_norm(1,i) = norm(V(1:3,i));
end
Vang_norm = zeros(1,n+1);
for i=1:n+1
    Vang_norm(1,i) = norm(V(4:6,i));
end
Plin_norm = zeros(1,n+1);
for i=1:n+1
    Plin_norm(1,i) = norm(P(1:3,i));
end
Pang_norm = zeros(1,n+1);
for i=1:n+1
    Pang_norm(1,i) = norm(P(4:6,i));
end
    
%figure
%a = plot3(V(1,:),V(2,:),V(3,:)),xlabel('x'),ylabel('y'),zlabel('z'),title('Linear Velocity');
%figure
%b = plot3(V(4,:),V(5,:),V(6,:)),xlabel('x'),ylabel('y'),zlabel('z'),title('Angular Velocity');
figure % Linear Velocity
plot(T,V(1,:),T,V(2,:),T,V(3,:));
xlabel('Time'),ylabel('Velocity'),title('Linear Velocity'),legend('X','Y','Z');
figure % Linear Velocity Norm
plot(T,Vlin_norm);
xlabel('Time'),ylabel('Velocity (Normed)'),title('Linear Velocity Norm');
% figure % Linear Momentum Norm
% plot(T,Plin_norm);
% xlabel('Time'),ylabel('Momentum (Normed)'),title('Linear Momentum Norm');
figure % Angular Velocity
plot(T,V(4,:),T,V(5,:),T,V(6,:));
xlabel('Time'),ylabel('Velocity'),title('Angular Velocity'),legend('X','Y','Z');
figure % Angular Velocity Norm
plot(T,Vang_norm);
xlabel('Time'),ylabel('Velocity (Normed)'),title('Angular Velocity Norm');
% figure % Angular Momentum Norm
% plot(T,Vang_norm);
% xlabel('Time'),ylabel('Momentum (Normed)'),title('Angular Momentum Norm');
