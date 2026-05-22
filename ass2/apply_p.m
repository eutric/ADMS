function F = apply_p(pG, elements, nodes, gamma, l, idb, incidenze)
lambda = @(g) [
    cos(g), sin(g), 0;
    -sin(g), cos(g), 0;
    0,0,1
];
Lambda = @(g) [
    lambda(g), zeros(3);
    zeros(3), lambda(g)
];
for ii=1:length(elements)
    gammaii = gamma(elements(ii)); % gamma of ii element
    lii = l(elements(ii));
    pL = lambda(gammaii)*pG;
    py = pL(2);
    F_kL = [
        0
        -py*lii/2;
        -py^2*lii/12;
        0;
        -py*lii/2;
        py^2*lii/12;
    ];% This is the vector of element ii, 6x1 e
    F_kG(:, ii) = Lambda*F_kL;
    Eii = zeros()
end
F = sum(F_kG);
end