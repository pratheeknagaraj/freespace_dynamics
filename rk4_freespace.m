function [R,P,V,T] = rk4_freespace( a,b,p0,v0,f_ext,g_ext,M,n )
%RK4_FREESPACE Modified RK4 for freespace equations
%   
%Modification to the RK4 routine that includes a subroutine
%for computing the momentum and velocity vectors
%
%Input   a - left endpoint
%        b - right endpoint
%        p0 - initial augmented (linear & angular) momentum vector
%        v0 - initial augmented (linear & angular) velocity vector
%        M - system inertia matrix
%        n - number of steps
%Ouput   p - final augmented (linear & angular) momentum vector
%        v - final augmented (linear & angular) velocity vector

h = (b-a)/n;        % Step Size
T = zeros(1,n+1);   % Time Vector
P = zeros(6,n+1);   % Momentum Vector Array
V = zeros(6,n+1);   % Velocity Vector Array
T = a:h:b;          % Populate Time Vector with intermediate points
P(:,1) = p0;        % Set initial momentum condition
V(:,1) = v0;        % Set initial velocity condition
M_inv = inv(M);     % Inverse System Intertia Matrix

for j = 1:n % Define loop to iterate over steps
    % Momementum Differential Equation
    Pdot = @(t,P) [-skew_sym(V(4:6,j)) zeros(3); -skew_sym(V(1:3,j)) -skew_sym(V(4:6,j))] * P + [f_ext(t); g_ext(t)];
        
    k1 = h * Pdot(T(j),P(:,j));           % Stage 1
    k2 = h * Pdot(T(j)+h/2,P(:,j)+k1/2);  % Stage 2
    k3 = h * Pdot(T(j)+h/2,P(:,j)+k2/2);  % Stage 3
    k4 = h * Pdot(T(j)+h,P(:,j)+k3);      % Stage 4
    
    P(:,j+1) = P(:,j) + (k1+2*k2+2*k3+k4)/6;
    V(:,j+1) = M_inv*P(:,j+1);
end

R = [T' P' V'];
