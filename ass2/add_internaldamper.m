function C = add_internaldamper(node_1, node_2, C, c, G, idb)
% Adds in global matrix K a spring from node_1 to node_2 of stiffness k
    x12 = [1,0,0,-1,0,0];
    Ck_local = c*(x12'*x12);
    lambda = @(G) [
        cos(G), sin(G), 0;
        -sin(G), cos(G), 0;
        0, 0, 1
    ];
    Lambda = [
        lambda(G), zeros(3);
        zeros(3), lambda(G)
    ];
    Ck_global = Lambda'*Ck_local*Lambda;
    i_dofs_k = [idb(node_1,:), idb(node_2,:)];
    C(i_dofs_k, i_dofs_k) = C(i_dofs_k, i_dofs_k) + Ck_global;
end