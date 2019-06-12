% CasADi SX variables not compatible with Matlab mean function. This
% implementation fixes the issue
function mean = casadi_mean(numbers)
    mean = sum1(numbers)/size(numbers,1);
end