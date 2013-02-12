function [R,P,V,T] = grk4_freespace( a,b,p0,v0,f_ext,g_ext,M,n )
%GRK4_FREESPACE Modified GRK4 for freespace equations
%   
%Modification to the Implicit Gauss-Legendre Order 4 Method
%routine that includes a subroutine for computing the momentum and velocity
%vectors. Both an algebraic/analytical and numerical subroutine is provided.
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

% GRK4 Butcher Tableau Values
a11 = 0.25;
a12 = 0.25 - sqrt(3)/6;
a21 = 0.25 + sqrt(3)/6;
a22 = 0.25;

b1 = 0.5;
b2 = 0.5;

c1 = 0.5 - sqrt(3)/6;
c2 = 0.5 + sqrt(3)/6;

for j = 1:n % Define loop to iterate over steps
    %Analytical/Algebraic Approach
%     A = [ -skew_sym(V(4:6,j)) zeros(3); -skew_sym(V(1:3,j)) -skew_sym(V(4:6,j)) ];
%     F1 = [ f_ext(T(j)+h*c1); g_ext(T(j)+h*c1) ];
%     F2 = [ f_ext(T(j)+h*c1); g_ext(T(j)+h*c2) ];
%     
%     k1 = inv( eye(6)-h*a11*A-h*h*a12*a21*A*inv(eye(6)-h*a22*A)*A )*(A*(P(:,j)+h*a12*inv(eye(6)-h*a22*A)*(A*P(:,j)+F2))+F1);
%     k2 = inv( eye(6)-h*a22*A-h*h*a21*a12*A*inv(eye(6)-h*a11*A)*A )*(A*(P(:,j)+h*a21*inv(eye(6)-h*a11*A)*(A*P(:,j)+F1))+F2);
%     
%     P(:,j+1) = P(:,j) + h*(b1*k1+b2*k2);
%     V(:,j+1) = M_inv*P(:,j+1);
    
    %Numerical Approach
    A = [ -skew_sym(V(4:6,j)) zeros(3); -skew_sym(V(1:3,j)) -skew_sym(V(4:6,j)) ];
    F1 = [ f_ext(T(j)+h*c1); g_ext(T(j)+h*c1) ];
    F2 = [ f_ext(T(j)+h*c2); g_ext(T(j)+h*c2) ];
    F = [F1; F2];
    Ap = A*P(:,j);
    Coeff = h*[a11*A a12*A; a21*A a22*A] - eye(12);
    
    func = @(K) [Ap; Ap] + Coeff*K + F;
    D_func = @(K) Coeff;
    K = newton_raphson( func, D_func, [P(:,j); P(:,j)] );
    k1 = K(1:6);
    k2 = K(7:12);
    
    P(:,j+1) = P(:,j) + h*(b1*k1 + b2*k2);
    V(:,j+1) = M_inv*P(:,j+1);
end

R = [T' P' V'];
