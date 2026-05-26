function F_G = F_t(F, x_t, t, n, ndof, l, gamma, xy, incidenze)    % 1 - Understand where i am
% This only works for vertical loads;
    xii = x_t(t);
    n_el = whereami([xii, 0], xy); % this time it's all in y=0
    l = l(n_el);
    a = xii-xy(n_el,1); % Again, only because everything is on y=0
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
        0;
        (F_L(2)*a*b*(b-a))/l^3 + F_L(2)*b/l;
        F_L(2)*a*b^2/l^2;
        0;
        (F_L(2)*a*b*(a-b))/l^3 + F_L(2)*a/l;
        F_L(2)*a^2*b/l^2
    ];
    
    Fk_global = Lambda(gamma(n_el))*Fk_local;
    E = zeros(length(Fk_local),n);
    E(:,incidenze(n_el,:)) = eye(length(Fk_local));
    F_G = E'*Fk_global;
    F_G = F_G(1:ndof);
    % Q_t = Phi'*F_G;
end
