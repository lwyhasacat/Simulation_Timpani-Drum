%  helper function for cotangent laplacian matrix using pre-generated mesh
%  [cotan, nodes]=getAngle(p, t, v)
%
%  p: node positions (x, y, z)
%  t: triangle indices or nodes number (#, #, #)
%  v: vertex shared by a list of triangles
% 
%  cotan: a list of cotangent of triangles in the list
%  nodes: indices of triangles in the list corresponding to cotangent list

function [cotan, nodes]=getAngle(p, t, v)

% finding triangles that share the same vertex
% use row number later for finding coordinates given by p
[row, column] = find(t == v);

% create a matrix to hold cotangent and nodes of triangles 
cotan = zeros(size(row,1), 2);
nodes = zeros(size(row,1), 3);

% looping through all the triangles (# given by the # of rows found)
for i = 1:size(row, 1)
    % column(i) is the column number of the vertex index
    % set col1 and col2 to column number of the two other triangle indices
    if column(i) == 1
        col1 = 2;
        col2 = 3;
    elseif column(i) == 2
        col1 = 1;
        col2 = 3;
    else
        col1 = 1;
        col2 = 2;
    end
    
    % getting coordinates for each triangle; each with x, y, z coordinates
    a_x = p(v, 1);
    a_y = p(v, 2);
    a_z = p(v, 3);
    b_x(i) = p(t(row(i), col1), 1);
    b_y(i) = p(t(row(i), col1), 2);
    b_z(i) = p(t(row(i), col1), 3);
    c_x(i) = p(t(row(i), col2), 1);
    c_y(i) = p(t(row(i), col2), 2);
    c_z(i) = p(t(row(i), col2), 3);
    
    % finding the vectors given by the triangles
    a(i,:) = [b_x(i)-c_x(i) b_y(i)-c_y(i) b_z(i)-c_z(i)];
    b(i,:) = [a_x-c_x(i) a_y-c_y(i) a_z-c_z(i)];
    c(i,:) = [a_x-b_x(i) a_y-b_y(i) a_z-b_z(i)];
    
    % get cotangent matrix, # of triangles * 2, using standard formula
    cotan(i, 1) = dot(a(i,:), b(i,:)) / sqrt(dot(a(i,:), a(i,:))...
                  * dot(b(i,:), b(i,:)) - dot(a(i,:), b(i,:))^2);
    cotan(i, 2) = dot(a(i,:), c(i,:)) / sqrt(dot(a(i,:), a(i,:))...
                  * dot(c(i,:), c(i,:)) - dot(a(i,:), c(i,:))^2);
end

% get corresponding nodes
for j = 1:size(row, 1)
    nodes(j, :) = t(row(j), :);
end
    
end