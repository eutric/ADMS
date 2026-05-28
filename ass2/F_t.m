function F_G = F_t(F, x_t, t, n, ndof, l, gamma, incidenze)   
% Time law data
x = x_t(t);
cum = [0, cumsum(l(23:30))];
n_el = 22 + find(x >= cum(1:end-1) & x < cum(2:end), 1, 'first');
l = l(n_el);
a = x - cumsum(n_el-22);
b = l-a;

lambda = @(g) [ 
    cos(g), -sin(g), 0;
    sin(g),  cos(g), 0;
    0,       0,      1;
];

Lambda = @(g) [ 
    lambda(g), zeros(3);
    zeros(3),  lambda(g)
];

F_L = lambda(gamma(n_el))'*F;
Fk_local = [
    F_L(1)*b/l;
    (F_L(2)*a*b*(b-a))/l^3 + F_L(2)*b/l;
    F_L(2)*a*b^2/l^2;
    F_L(1)*a/l;
    (F_L(2)*a*b*(a-b))/l^3 + F_L(2)*a/l;
    -F_L(2)*a^2*b/l^2
];

Fk_global = Lambda(gamma(n_el))*Fk_local;
E = zeros(length(Fk_local),n);
E(:,incidenze(n_el,:)) = eye(length(Fk_local));
F_G = E'*Fk_global;
F_G = F_G(1:ndof);
end
