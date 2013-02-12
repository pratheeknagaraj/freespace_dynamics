function [R,P,V,T] = impliciteuler_freespace( a,b,p0,v0,f_ext,g_ext,M,n )
%IMPLICITEULER_FREESPACE Modified Implicit (Backwards) Euler Method for freespace equations
%   
%Modification to the Implicit (Backwards) Euler routine that includes a subroutine
%for computing the momentum and velocity vectors. 3 subroutines are
%provided in the loop, with 2 algebriac approaches and 1 numerical approach
%by implementing the Newton-Raphson method.
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
    % ======= Algebraic Approach Start =======
    % Solve Algebraic Equation and compute new vectors
    % NOTE: Algebraic Manipulation works, although steps must be large to see
    % accuracy
    % -- Uncomment Region 1 Start --
%     p_new = inv( eye(3) + h * skew_sym(V(4:6,j)) )*( P(1:3,j) + h * f_ext(T(j)+h));
%     h_new = inv( eye(3) + h * skew_sym(V(4:6,j)) )*( P(4:6,j)-h*(skew_sym(V(1:3,j))*p_new + g_ext(T(j)+h)) );
%     
%     P(:,j+1) = [ p_new; h_new ];
%     V(:,j+1) = M_inv*P(:,j+1);
    % --- Uncomment Region 1 End ---
    %
    % -- Uncomment Region 2 Start --
%     P(:,j+1) = inv( eye(6)-h*[ -skew_sym(V(4:6,j)) zeros(3); -skew_sym(V(1:3,j)) -skew_sym(V(4:6,j))] )*( P(:,j) + h*[ f_ext(T(j)+h); g_ext(T(j)+h) ] );
%     V(:,j+1) = M_inv*P(:,j+1);
    % --- Uncomment Region 2 End ---
    % ======== Algebraic Approach End ========
    
    % ======= Numerical Approach Start =======
    % Solve Numerical Equation and compute new vectors
    coeff = ( eye(6)-h*[ -skew_sym(V(4:6,j)) zeros(3); -skew_sym(V(1:3,j)) -skew_sym(V(4:6,j))] );
    func = @(P_new) coeff*P_new - P(:,j) - h*[ f_ext(T(j)+h); g_ext(T(j)+h) ];
    D_func = @(P_new) coeff;
    P(:,j+1) = newton_raphson( func, D_func, P(:,j) );
    V(:,j+1) = M_inv*P(:,j+1);
    
end

R = [T' P' V'];
