function [R,P,V,T,steps] = dopri54_freespace( a,b,p0,v0,f_ext,g_ext,M,tol )
%GRK4_FREESPACE Modified GRK4 for freespace equations
%   
%Modification to the Implicit Gauss-Legendre Order 4 Method
%routine that includes a subroutine for computing the momentum and velocity
%vectors. Only a algebraic/analytical subroutine is provided.
%
%Input   a - left endpoint
%        b - right endpoint
%        p0 - initial augmented (linear & angular) momentum vector
%        v0 - initial augmented (linear & angular) velocity vector
%        M - system inertia matrix
%        tol - error tolerance per step
%Ouput   p - final augmented (linear & angular) momentum vector
%        v - final augmented (linear & angular) velocity vector

h = tol^(1/5)/4;    % Step Size
T = zeros(1,1);     % Time Vector
P = zeros(6,1);     % Momentum Vector Array
V = zeros(6,1);     % Velocity Vector Array
T = zeros(1);       % Populate Time Vector with intermediate points
P(:,1) = p0;        % Set initial momentum condition
V(:,1) = v0;        % Set initial velocity condition
T(1) = a;           % Set initial time
M_inv = inv(M);     % Inverse System Intertia Matrix

% GRK4 Butcher Tableau Values
a2 = [1/5];
a3 = [3/40 9/40];
a4 = [44/45 -56/15 32/9];
a5 = [19372/6561 -25360/2187 64448/6561 -212/729];
a6 = [9017/3168 -355/33 46732/5247 49/176 -5103/18656];
a7 = [35/384 0 500/1113 125/192 -2187/6784 11/84];

b1 = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];
b2 = [35/384 0 500/1113 125/192 -2187/6784 11/84 0];

c = [0 1/5 3/10 4/5 8/9 1 1];

steps = 0;
nrej = 0;
j = steps+1;
while T(j) < b
    if T(j)+h > b
        h = b-T(j);
    end
    Pdot = @(t,P) [-skew_sym(V(4:6,j)) zeros(3); -skew_sym(V(1:3,j)) -skew_sym(V(4:6,j))] * P + [f_ext(t); g_ext(t)];
        
    k1 = h * Pdot(T(j)       ,P(:,j));                                                          % Stage 1
    k2 = h * Pdot(T(j)+h*c(2),P(:,j)+k1*a2(1));                                                 % Stage 2
    k3 = h * Pdot(T(j)+h*c(3),P(:,j)+k1*a3(1)+k2*a3(2));                                        % Stage 3
    k4 = h * Pdot(T(j)+h*c(4),P(:,j)+k1*a4(1)+k2*a4(2)+k3*a4(3));                               % Stage 4
    k5 = h * Pdot(T(j)+h*c(5),P(:,j)+k1*a5(1)+k2*a5(2)+k3*a5(3)+k4*a5(4));                      % Stage 5
    k6 = h * Pdot(T(j)+h*c(6),P(:,j)+k1*a6(1)+k2*a6(2)+k3*a6(3)+k4*a6(4)+k5*a6(5));             % Stage 6
    k7 = h * Pdot(T(j)+h*c(7),P(:,j)+k1*a7(1)+k2*a7(2)+k3*a7(3)+k4*a7(4)+k5*a7(5)+k6*a7(6));    % Stage 7
    
    Y5 = b1(1)*k1+b1(2)*k2+b1(3)*k3+b1(4)*k4+b1(5)*k5+b1(6)*k6+b1(7)*k7;
    Y4 = b2(1)*k1+b2(2)*k2+b2(3)*k3+b2(4)*k4+b2(5)*k5+b2(6)*k6+b2(7)*k7;
    err = norm(Y5-Y4);
    
    if err < tol
        T(j+1) = T(j) + h;
        P(:,j+1) = P(:,j) + Y5;
        V(:,j+1) = M_inv*P(:,j+1);
        steps = steps+1;
        j = j+1;
    else
        nrej = nrej+1;
    end
    
    delta = 0.84*(tol/err)^(1/5);
    if delta < 0.1
        h = h * 0.1;
    elseif delta > 4.0
        h = h * 4.0;
    else
        h = delta * h;
    end   
end

R = [T' P' V'];
