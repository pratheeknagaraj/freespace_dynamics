function [ skew ] = skew_sym( c )
%Skew_Sym Produces a skew symmetric matrix of input vector
%   Forms a skew-symmetric matrix with input vector
skew = [0 -1*c(3) c(2); c(3) 0 -1*c(1); -1*c(2) c(1) 0]; %#ok<NOPRT>

end
