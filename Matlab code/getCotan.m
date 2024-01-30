%  helper function for cotangent laplacian matrix using pre-generated mesh
%  [sumCotan, node_list_e]=getCotan(p, t, v)
%
%  p: node positions (x, y, z)
%  t: triangle indices or nodes number (#, #, #)
%  v: vertex shared by a list of triangles that we want to look at
% 
%  sumCotan: sum of cotangent of theta_1,2 in T1 and T2 sharing an edge
%  nodes_list_e: a list of indices of nodes adjacent to the vertex

function [sumCotan, node_list_e]=getCotan(p, t, v)

% get a list of the edges on the boundary
% all possible edges are given by combinations of column 1, 2, and 3
edges = [t(:,[1,2]);
        t(:,[1,3]);
        t(:,[2,3])];
% sort the edges
edges = sort(edges,2);
% ~ is supposedly a matrix C that contains only the unique rows from edges
% edges(:) = C(ic)
[~,ia,ic] = unique(edges,'rows');
% check how many times ic(i) appears, i as indexing var
% if ic only appears one time, this edge only appeared one time
% and therefore a boundary edge because only one triangle contains it
count = histc(ic,1:max(ic));
num = find(count==1);
% retrieving the boundary edges
be = edges(ia(num),:);

% finding triangles that share the same vertex
% use row number later for finding coordinates given by p
[row, ~] = find(t == v); % store row/column number where we find v
n = size(row,1); % number of triangles that are important here
node_list = zeros(n, 3); % a list for storing node indices
node_list_vec = zeros(n*3, 1); % a list later used for special case
node_list_e = zeros(n*3, 1); % a list for storing adjacent node indices

% get node indices of the triangles we care about here
for r = 1:n
    for c = 1:3
        node_list(r, c) = t((row(r)), c);
    end
end

% get all node indices adjacent to the original vertex v (excluding itself)
count = 1;
for r = 1:n
    for c = 1:3
        node_list_e(count) = t(row(r), c);
        node_list_vec(count) = t(row(r), c);
        count = count + 1;
    end
end
node_list_e = setdiff(unique(node_list_e), v);
node_list_e_dup = node_list_e; % later used after modifying node_list_e

% check if any edge given by v and node_list_e(j) is a boundary edge
% if so, delete node_list_e(j) from node_list_e
% because those edges are not shared by any of the two triangles
store_list = [];
s = size(node_list_e, 1);
% loop through my node_list_e to check whether any of the edge matches with
% my boundary edge list be that's found earlier
for j = 1:s
    v_ej = sort([v, node_list_e(j)], 2);
    checker = ismember(v_ej, be, 'rows');
    if checker == 1
        % store the location where I find the boundary edges in node_list_e
        store_list(end+1) = j; 
    end
end
% deleting those boundary edges from my node_list_e
% since I'm deleting elements, the index of the boundary edges is changing 
% -k+1 takes care of this indexing problem 
for k = 1:size(store_list, 2)
    node_list_e(store_list(1, k)-k+1) = [];
end

% this if section of my code deals with the case where there are only two
% triangles instead of a ring (case 2 in the paper)
% the logic of this section is still the same as the one dealing with a
% ring of triangles
if numel(node_list_e) == 1
    % edge_nodes: the node indices that forms the edge shared by T1, T2
    edge_nodes = [v node_list_e(1)];
    % tri_nodes: the the other two node indices in T1, T2
    tri_nodes = setdiff(node_list_e_dup, edge_nodes);
    
    % find the corresponding coordinates for nodes in T1, T2
    v_x = p(edge_nodes(1), 1);
    v_y = p(edge_nodes(1), 2);
    v_z = p(edge_nodes(1), 3);
    
    vadj_x = p(edge_nodes(2), 1);
    vadj_y = p(edge_nodes(2), 2);
    vadj_z = p(edge_nodes(2), 3);
    
    vt1_x = p(tri_nodes(1), 1);
    vt1_y = p(tri_nodes(1), 2);
    vt1_z = p(tri_nodes(1), 3);
    
    vt2_x = p(tri_nodes(2), 1);
    vt2_y = p(tri_nodes(2), 2);
    vt2_z = p(tri_nodes(2), 3);
    
    % finding the vectors given by the triangles
    a = [v_x-vt1_x v_y-vt1_y v_z-vt1_z];
    b = [vadj_x-vt1_x vadj_y-vt1_y vadj_z-vt1_z];
    c = [v_x-vt2_x v_y-vt2_y v_z-vt2_z];
    d = [vadj_x-vt2_x vadj_y-vt2_y vadj_z-vt2_z];
    
    % get cotangent matrix, using standard formula
    cotan = [dot(b, a) / sqrt(dot(b, b) * dot(a, a) - dot(b, a)^2) ...
                dot(d, c) / sqrt(dot(d, d) * dot(c, c) - dot(d, c)^2)];
            
    % sumCotan
    sumCotan = sum(cotan, 2);
    return
end

% create matrices to hold cotan, nodes in T1 and T2, and sum of cotan
cotan = zeros(n, 2);
sumCotan = zeros(n, 1);

% create matrices to hold vectors later used for calculating cotangent
a = zeros(size(node_list_e, 1), 3);
b = zeros(size(node_list_e, 1), 3);
c = zeros(size(node_list_e, 1), 3);
d = zeros(size(node_list_e, 1), 3);

% loop through all the node indices adjacent to v
for i = 1:size(node_list_e, 1)
    % for finding T1, T2 that share the same edge (v and an adjacent node)
    % row_T and column_T give the row/col number of the adj node
    [row_T, column_T] = find(node_list == node_list_e(i));
    
    % co(1), co(2) are the col number of v corresponding to adj node
    co = [];
    co(1) = find(node_list(row_T(1), :) == v);
    co(2) = find(node_list(row_T(2), :) == v);
    
    % get the last unknown node v_t1 in T1, v_t2 T2
    % col(1) is the column # for v_t1 in T1
    % col(2) is the column # for v_t2 in T2
    
    % Explaining the logic:
    % since we have 2 triangles T1, T2 sharing edge e, we have 4 vertices
    % v, vi depending on our current loop <- 2 vertices shared by T1, T2,
    % and v_t1, v_t2. Here we are trying to find v_t1, v_t2. We already 
    % know the row number of t that tells us the 3 node indices in triangle
    % T1 and same for T2, and we've already found the column number for v, 
    % vi in T1 and T2. Since the matrix t has only 3 columns, we can figure
    % out the column number of v_t1 and v_t2 by taking a set difference
    % between all possible column numbers 1, 2, 3 and the known column
    % number column_T(t_i) and co(t_i); t_i is a looping index for storing
    % column number of v_t1 and column number of v_t2 into the list col
    col = [];
    for t_i = 1:2
        col_list = setdiff([1 2 3], column_T(t_i));
        col(t_i) = sum(setdiff(col_list, co(t_i)));
    end
    % summary for convenience:
    %       column_T,       co,       col
    %           |           |          |
    %       adj node,       v,        v_t1
    
    % find the corresponding coordinates for nodes in T1, T2
    v_x = p(v, 1);
    v_y = p(v, 2);
    v_z = p(v, 3);
    
    vadj_x(i) = p(node_list(row_T(1), column_T(1)), 1);
    vadj_y(i) = p(node_list(row_T(1), column_T(1)), 2);
    vadj_z(i) = p(node_list(row_T(1), column_T(1)), 3);
    
    vt1_x(i) = p(node_list(row_T(1), col(1)), 1);
    vt1_y(i) = p(node_list(row_T(1), col(1)), 2);
    vt1_z(i) = p(node_list(row_T(1), col(1)), 3);
    
    vt2_x(i) = p(node_list(row_T(2), col(2)), 1);
    vt2_y(i) = p(node_list(row_T(2), col(2)), 2);
    vt2_z(i) = p(node_list(row_T(2), col(2)), 3);
    
    % finding the vectors given by the triangles
    a(i,:) = [v_x-vt1_x(i) v_y-vt1_y(i) v_z-vt1_z(i)];
    b(i,:) = [vadj_x(i)-vt1_x(i) vadj_y(i)-vt1_y(i) vadj_z(i)-vt1_z(i)];
    c(i,:) = [v_x-vt2_x(i) v_y-vt2_y(i) v_z-vt2_z(i)];
    d(i,:) = [vadj_x(i)-vt2_x(i) vadj_y(i)-vt2_y(i) vadj_z(i)-vt2_z(i)];
    
    % get cotangent matrix, using standard formula
    cotan(i, 1) = dot(b(i,:), a(i,:)) / sqrt(dot(b(i,:), b(i,:))...
                  * dot(a(i,:), a(i,:)) - dot(b(i,:), a(i,:))^2);
    cotan(i, 2) = dot(d(i,:), c(i,:)) / sqrt(dot(d(i,:), d(i,:))...
                  * dot(c(i,:), c(i,:)) - dot(d(i,:), c(i,:))^2);
    
    % sumCotan
    sumCotan = sum(cotan, 2);
end

end
