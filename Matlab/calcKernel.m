function g = calcKernel(Z, b) 

% Function to evaluate the assortativity kernel g(z) with parameter b
%
% USAGE: g = calcKernel(Z, b) 
%
% INPUTS: Z - a scalar, vector, or matrix of inputs to the assortativity
% kerenl (these could represent the difference between the activity level
% quantiles (x-y) or the difference between the actual activity levels
% (vx-vy)
%         b - non-negative scalar parameter for the assortativity kernel
%         (larger b means stronger assortativity and b=0 should reduce to
%         proportionate mixing)
%
% OUTPUTS: g - kernel evaluated at each element in Z


g = exp(-b*Z.^2);

