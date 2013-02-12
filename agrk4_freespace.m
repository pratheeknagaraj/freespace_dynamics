function [R,P,V,T,steps] = agrk4_freespace( a,b,p0,v0,f_ext,g_ext,M,tol )
%AGRK4_FREESPACE Modified Adaptive GRK4 for freespace equations
%   
%Modification to the Adapative Implicit Gauss-Legendre Order 4 Method
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

h = tol^(1/4);    % Step Size
T = zeros(1,1);     % Time Vector
P = zeros(6,1);     % Momentum Vector Array
V = zeros(6,1);     % Velocity Vector Array
T = zeros(1);       % Populate Time Vector with intermediate points
P(:,1) = p0;        % Set initial momentum condition
V(:,1) = v0;        % Set initial velocity condition
T(1) = a;           % Set initial time
M_inv = inv(M);     % Inverse System Intertia Matrix

% GRK4 Butcher Tableau Values
a11 = 0.25;
a12 = 0.25 - sqrt(3)/6;
a21 = 0.25 + sqrt(3)/6;
a22 = 0.25;

b1 = 0.5;
b2 = 0.5;

c1 = 0.5 - sqrt(3)/6;
c2 = 0.5 + sqrt(3)/6;

steps = 0;
nrej = 0;
j = steps+1;
while T(j) < b
    if T(j)+h > b
        h = b-T(j);
    end
    % ############# Double GRK4 Approach ################# %
    % -------- Analytical/Algebraic Approach Start ------- %
    A = [ -skew_sym(V(4:6,j)) zeros(3); -skew_sym(V(1:3,j)) -skew_sym(V(4:6,j)) ];
    F1 = [ f_ext(T(j)+h*c1); g_ext(T(j)+h*c1) ];
    F2 = [ f_ext(T(j)+h*c2); g_ext(T(j)+h*c2) ];

    k1 = inv( eye(6)-h*a11*A-h*h*a12*a21*A*inv(eye(6)-h*a22*A)*A )*(A*(P(:,j)+h*a12*inv(eye(6)-h*a22*A)*(A*P(:,j)+F2))+F1);
    k2 = inv( eye(6)-h*a22*A-h*h*a21*a12*A*inv(eye(6)-h*a11*A)*A )*(A*(P(:,j)+h*a21*inv(eye(6)-h*a11*A)*(A*P(:,j)+F1))+F2);

    Y4 = P(:,j) + h*(b1*k1+b2*k2);
    
    Vi = V(:,j);
    Pi = P(:,j);
    h = h/2;
    for i = 1:2
        A = [ -skew_sym(Vi(4:6,i)) zeros(3); -skew_sym(Vi(1:3,i)) -skew_sym(Vi(4:6,i)) ];
        F1 = [ f_ext(T(j)+i*h*c1); g_ext(T(j)+i*h*c1) ];
        F2 = [ f_ext(T(j)+i*h*c2); g_ext(T(j)+i*h*c2) ];

        k1 = inv( eye(6)-h*a11*A-h*h*a12*a21*A*inv(eye(6)-h*a22*A)*A )*(A*(Pi(:,i)+h*a12*inv(eye(6)-h*a22*A)*(A*Pi(:,i)+F2))+F1);
        k2 = inv( eye(6)-h*a22*A-h*h*a21*a12*A*inv(eye(6)-h*a11*A)*A )*(A*(Pi(:,i)+h*a21*inv(eye(6)-h*a11*A)*(A*Pi(:,i)+F1))+F2);

        Pi(:,i+1) = Pi(:,i) + h*(b1*k1+b2*k2);
        Vi(:,i+1) = M_inv*Pi(:,i+1);
    end
    h = h*2;
    
    err = norm(Pi(:,3)-Y4);
    
    if err < tol
        T(j+1) = T(j) + h;
        P(:,j+1) = Pi(:,3);
        V(:,j+1) = Vi(:,3);
        steps = steps+1;
        j = j+1;
    else
        nrej = nrej+1;
    end
    
    delta = 0.9*(tol/err)^(1/4);
    if delta < 0.1
        h = h * 0.1;
    elseif delta > 4.0
        h = h * 4.0;
    else
        h = delta * h;
    end
    
    % -------- Analytical/Algebraic Approach End -------- %
    
    % Numerical Approach
%     A = [ -skew_sym(V(4:6,j)) zeros(3); -skew_sym(V(1:3,j)) -skew_sym(V(4:6,j)) ];
%     F1 = [ f_ext(T(j)+h*c1); g_ext(T(j)+h*c1) ];
%     F2 = [ f_ext(T(j)+h*c1); g_ext(T(j)+h*c2) ];
%     F = [F1; F2];
%     Ap = A*P(:,j);
%     Coeff = h*[a11*A a12*A; a21*A a22*A] - eye(12);
%     
%     func = @(K) [Ap; Ap] + Coeff*K + F;
%     D_func = @(K) Coeff;
%     K = newton_raphson( func, D_func, [P(:,j); P(:,j)] );
%     k1 = K(1:6);
%     k2 = K(7:12);
%     
%     P(:,j+1) = P(:,j) + h*(b1*k1 + b2*k2);
%     V(:,j+1) = M_inv*P(:,j+1);
    % -- End Double GRK4 Approach -- %
end

R = [T' P' V'];
