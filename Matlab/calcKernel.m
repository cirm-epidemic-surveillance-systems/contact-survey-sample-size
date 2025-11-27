function g = calcKernel(Z, b) 

% Evaluate the kernel g(z) for the assortativity model with parameter b

g = exp(-b*Z.^2);

