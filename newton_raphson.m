function x0 = newton_raphson( f, fprime, start )
% Newton-Raphson method to iteratively compute zeros of function

x0 = start;                                 % Initial Guess
zero = ones(size(start,1),1);               % Vector to hold final answer
steps = 0;                                  % Keep track of iterations
threshold = ones(size(start,1),1)*1e-12;    % Error threshold

while (abs(zero) > threshold)               % Iteration loop for Newton-Raphson method
    zero = inv(fprime(x0))*(f(x0));
    x0 = x0 - zero;
    steps = steps + 1;
    if (steps > 100)                        % Break if not converging
        disp('Error in convergence')
        x0 = 0;
        break
    end
end
