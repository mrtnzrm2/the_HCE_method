function partition = fastCutTree(H, n_clusters)
    % Custom cut_tree implementation for hierarchical clustering
    N = size(H,1) + 1;
    T = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    for i = 1:N
        T(i) = i;
    end

    K = N;
    i = 1;

    while true
        if K <= n_clusters
            break;
        end

        nx = H(i, 1);
        ny = H(i, 2);

        T(N+i) = [T(nx), T(ny)];

        remove(T, nx);
        remove(T, ny);

        i = i + 1;
        K = K - 1;
    end

    partition = zeros(N,1);
    keysT = T.keys;
    for k = 1:length(keysT)
        idx = keysT{k};
        nodes = T(idx);
        partition(nodes) = idx;
    end
end