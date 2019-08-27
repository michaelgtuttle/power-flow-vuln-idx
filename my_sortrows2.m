function [sorted_matrix] = my_sortrows(matrix, col)

[m,n] = size(matrix);
sorted_matrix = zeros(m, n);
cur_line = 1;

while size(matrix, 1) ~= 0 
    max_idx = 1;
    max_val = matrix(max_idx, col);    
    for i = 1:size(matrix, 1)
        if matrix(i, col) > max_val
            max_val = matrix(i, col);
            max_idx = i;
        end
    end
    sorted_matrix(cur_line, :) = matrix(max_idx, :);
    matrix(max_idx, :) = [];
    cur_line = cur_line + 1;
end
    


