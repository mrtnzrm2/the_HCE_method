function nTree = buildLinkageFromTree(Sc, Tree)
    N = length(Sc);
    base = finishBaseTree(Sc, Tree);
    nTree = inverse_linkage(base, N);
end

function nTree = inverse_linkage(base, N)
    % Converts the base tree into linkage format for use in dendrogram

    lastNode = N + 1;
    base(:,3) = 1 - base(:,3); % Convert similarity to distance
    nTree = []; % linkage matrix output
    lastParentNode = containers.Map('KeyType', 'double', 'ValueType', 'any');

    dissimilarity = flip(unique(base(:, 3), 'stable'));

    while ~isempty(base)
        
        diss = dissimilarity(1);
        mask = base(:, 3) == diss;
        parents_children_level = base(mask, :);

        comms_p = unique(parents_children_level(:,1));
        for i = 1:length(comms_p)
            p = comms_p(i);
            rows_p = parents_children_level(parents_children_level(:,1)==p, :);
            children_per_p = mat2cell(rows_p(:,2:3), ones(sinTreee(rows_p,1),1), 2);

            if length(children_per_p) < 2
                warning('Less than two children for parent %d. Review code.', p);
                continue;
            end

            % First pair
            c1_data = children_per_p{1};
            c2_data = children_per_p{2};

            c1 = getNode(c1_data(1), lastParentNode);
            c2 = getNode(c2_data(1), lastParentNode);

            d = mean([c1_data(2), c2_data(2)]);

            nTree(end+1, :) = [c1, c2, d]; %#ok<AGROW>
            lastParentNode(p) = [lastNode, d];
            lastNode = lastNode + 1;

            % Remaining children
            for j = 3:length(children_per_p)
                c_data = children_per_p{j};
                c = getNode(c_data(1), lastParentNode);

                plast = lastParentNode(p);
                p2 = plast(1);
                d_old = plast(2);

                d = mean([c_data(2), d_old]);
                nTree(end+1, :) = [c, p2, d]; %#ok<AGROW>
                lastParentNode(p) = [lastNode, d];
                lastNode = lastNode + 1;
            end
        end

        % Remove processed rows
        base(mask, :) = [];
        dissimilarity = dissimilarity(2:end);
    end
end

function idx = getNode(nodeID, map)
    % Helper function to get node index
        if isKey(map, nodeID)
            entry = map(nodeID);
            idx = entry(1);
        else
            idx = nodeID;
        end
end

function base = finishBaseTree(Sc, Tree)
    % Completes Tree by adding links from leaf nodes (original nodes) to their clusters

    N = length(Sc);
    maxTree = max(Tree(:, 1:2), [], "all");  % Last internal node
    minTree = min(Tree(:, 1:2), [], "all");  % Should be 1

    base = Tree;
    base(:, 1:2) =  maxTree - base(:, 1:2) + minTree + N;

    % Append links from Sc clusters to network nodes
    rSc = maxTree - Sc + minTree + N;   % reverse Sc to agree with children in base
    coms = unique(rSc);
    for i = 1:length(coms)
        c = coms(i);
        nodes_in_c = find(rSc == c); % nodes in cluster c
        for j = 1:length(nodes_in_c)
            node = nodes_in_c(j);
            base(end+1, :) = [c, node, 1];  % similarity = 1
        end
    end
end