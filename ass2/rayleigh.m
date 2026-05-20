function x = rayleigh(om, xi)
% Function that returns alpha and beta rayleigh damping coefficients, from
% a known vector of damping of the firsts length(xi) modes
% the vector of corresponding natural freq is needed

    AA = zeros(length(xi), 2);
    for ii=1:length(xi)
        AA(ii,:) = [1/2/om(ii), om(ii)/2];
    end
    A = AA'*AA;
    b = AA'*xi;
    
    x = A\b;
end