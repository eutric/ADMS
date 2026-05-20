function K = add_spring(node_1, node_2, K, k, G, idb)
% Adds in global matrix K a spring from node_1 to node_2 of stiffness k
    x12 = [1,0,0,-1,0,0];
    Kk_local = k*(x12'*x12);
    lambda = @(G) [
        cos(G), sin(G), 0;
        -sin(G), cos(G), 0;
        0, 0, 1
    ];
    Lambda = [
        lambda(G), zeros(3);
        zeros(3), lambda(G)
    ];
    Kk_global = Lambda'*Kk_local*Lambda;
    i_dofs_k = [idb(node_1,:), idb(node_2,:)];
    K(i_dofs_k, i_dofs_k) = K(i_dofs_k, i_dofs_k) + Kk_global;
end