function [labels, K] = findHCELevel(H, nK, rn)
% Find best partition from linkage matrix using HCE
% H is (N-1)x3
% rn is the number of renormalization rounds

    if nargin < 2
        nK = [];
    end
    if nargin < 3
        rn = 0;
    end

    N = size(H, 1) + 1;

    if isempty(nK)
        hce = HCE(H, N);
        [K, ~] = getBestHCELevel(hce);
        for i = 1:rn
            hce = rHCE(H, N, K);
            [K, ~] = getBestHCELevel(hce);
        end
    else
        K = nK;
    end

    labels = fastCutTree(H, K);
end

function hce = HCE(H, N)
    
    % Compute Hierarchical Clustering Entropy across all levels
    % H is the (N-1)x3 linkage matrix
    % N is the number of original observations

    T = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    for i = 1:N
        T(i) = struct('size', 0);
    end

    hce = containers.Map('KeyType', 'int32', 'ValueType', 'double');

    for i = 1:(N-1)
        nx = int64(H(i, 1));
        ny = int64(H(i, 2));
        h = H(i, 3);

        T(N + i) = struct('size', T(nx).size + T(ny).size + 1);
        
        remove(T, nx);
        remove(T, ny);

        s = 0;
        keysT = T.keys;
        for j = 1:length(keysT)
            key = keysT{j};
            sz = T(key).size;
            if sz == 0
                continue;
            end
            p = double(sz) / double(i);
            s = s - p * log(p);
        end

        if i < N-1
            if h ~= H(i+1, 3)
                hce(N - i) = s * double(i) / double(N - 1);
            end
        else
            hce(N - i) = s * double(i) / double(N - 1);
        end
    end
end

function hce = rHCE(H, N, Kscale)
    % Renormalized HCE calculation
    T = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    for i = 1:N
        T(i) = struct('size', 0);
    end

    hce = containers.Map('KeyType', 'int32', 'ValueType', 'double');

    % Pre-merge steps to renormalize
    for i = 1:(N - Kscale)
        nx = int64(H(i, 1));
        ny = int64(H(i, 2));

        T(N + i) = struct('size', 0);
        remove(T, nx);
        remove(T, ny);

        hce(N - i) = 0;
    end

    % Actual rHCE computation
    for j = 1:(Kscale - 1)
        i = N - Kscale + j;
        nx = int64(H(i, 1));
        ny = int64(H(i, 2));
        h = H(i, 3);

        T(N + i) = struct('size', T(nx).size + T(ny).size + 1);

        remove(T, nx);
        remove(T, ny);

        s = 0.;
        keysT = T.keys;
        for k = 1:length(keysT)
            key = keysT{k};
            sz = T(key).size;
            if sz == 0
                continue;
            end
            p = double(sz) / double(j);
            s = s - p * log(p);
        end

        if i < N-1
            if h ~= H(i+1, 3)
                hce(N - i) =  s * double(j) / double(Kscale - 1);
            end
        else
            hce(N - i) = s * double(j) / double(Kscale - 1);
        end
    end
end



function [realK, maxHCE] = getBestHCELevel(hceMap, f)
    % Select best K (number of communities) based on HCE values
    if nargin < 2
        f = [];
    end

    K = cell2mat(hceMap.keys);
    H2 = cell2mat(hceMap.values);
    [~, idx] = sort(K, 'ascend');
    K = K(idx);
    H2 = H2(idx);

    [~, maxIdx] = max(H2);
    realK = K(maxIdx);
    maxHCE = H2(maxIdx);

    if isa(f, 'function_handle')
        realK = f(realK);
    end
    realK = int64(realK);
end