function R=rk4(f,a,b,ya,M)
%Classical Runge-Kutta 4 Implementation
%
%Input   f - the function
%        a - left endpoint
%        b - right endpoint
%        ya - initial condition y(a)
%        M - number of steps
%Ouput   R = [X' Y'] - X is vector of absicissas
%                    - Y is vector of ordinates

h = (b-a)/M; % Step Size
X = zeros(1,M+1); % Absicissas Vector
Y = zeros(1,M+1); % Ordinate Vector
X = a:h:b; % Populate Absicissas Vector with intermediate points
Y(1) = ya; % Set initial condition
for j = 1:M % Define loop to iterate over steps
    k1 = h * f(X(j),Y(j)); % Stage 1
    % -- Insert Subprocess --
    k2 = h * f(X(j)+h/2,Y(j)+k1/2); % Stage 2
    % -- Insert Subprocess --
    k3 = h * f(X(j)+h/2,Y(j)+k2/2); % Stage 3
    % -- Insert Subprocess --
    k4 = h * f(X(j)+h,Y(j)+k3); % Stage 4
    % -- Insert Subprocess --
    Y(j+1) = Y(j)+(k1+2*k2+2*k3+k4)/6;
end

R = [X' Y'];

