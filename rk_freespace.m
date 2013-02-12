function [R,P,V,T] = rk_freespace( a,b,p0,v0,f_ext,g_ext,M,n,type )
%RK_FREESPACE Runge-Kutta Methods for freespace equations
%   
%Modification to the Runge-Kutta Method routine that includes a 
%subroutine for computing the momentum and velocity vectors. Only works for
%non-adaptive routines. Both explicit and implicit methods work. Numerical
%procedure using newton-raphson method to solve for stage constants (k).
%
%Input   a - left endpoint
%        b - right endpoint
%        p0 - initial augmented (linear & angular) momentum vector
%        v0 - initial augmented (linear & angular) velocity vector
%        M - system inertia matrix
%        n - number of steps
%        type - specific Runge-Kutta method to evaluate
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
RK_Method = RK_coefficients();
[a,b,c] = RK_Method.coef(type);
stages = length(c);

for j = 1:n % Define loop to iterate over steps   
    %Numerical Approach
    A = [ -skew_sym(V(4:6,j)) zeros(3); -skew_sym(V(1:3,j)) -skew_sym(V(4:6,j)) ];
    F = zeros(stages*6,1);
    for k = 1:stages
        F((k-1)*6+1:k*6,1) = [f_ext(T(j)+h*c(k)); g_ext(T(j)+h*c(k))];
    end
    Ap = ones(6*stages,1);
    for n = 1:stages
        Ap((n-1)*6+1:n*6,1) = A*P(:,j);
    end
    BigA = zeros(stages*6);
    for i1 = 1:stages
        for i2 = 1:stages
            BigA((i1-1)*6+1:i1*6,(i2-1)*6+1:i2*6) = a(i1,i2)*A;
        end
    end
    Coeff = h*BigA - eye(stages*6,stages*6);
    
    func = @(K) Ap + Coeff*K + F;
    D_func = @(K) Coeff;
    init = zeros(stages*6,1);
    for m = 1:stages
        init((m-1)*6+1:m*6,1) = P(:,j);
    end
    K = newton_raphson( func, D_func, init );
    weight = zeros(6,1);
    for l = 1:stages
        weight = weight + b(l)*K((l-1)*6+1:l*6);
    end
    
    P(:,j+1) = P(:,j) + h*weight;
    V(:,j+1) = M\P(:,j+1);
end

R = [T' P' V'];
