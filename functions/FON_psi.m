% Compute psi values for a vector alpha
% Inputs:
%   alpha - Fractional order vector
%   j     - Number of iterations
% Outputs:
%   psi   - Matrix (j × length(alpha)) containing psi values
%          for each alpha and iteration
function psi = FON_psi(alpha, j)

    psi(1,:) = 1;
    for i = 1:j
        psi(i+1,:) = psi(i,:) .* ((i-1-alpha)/i);
    end
    
end