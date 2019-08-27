function [sorted_matrix] = my_sortrows(matrix, col)

[m,n] = size(matrix);
sorted_matrix = zeros(m, n);
cur_line = 1;

while size(matrix, 1) ~= 0 
    min_idx = 1;
    min_val = matrix(min_idx, col);    
    for i = 1:size(matrix, 1)
        if matrix(i, col) < min_val
            min_val = matrix(i, col);
            min_idx = i;
        end
    end
    sorted_matrix(cur_line, :) = matrix(min_idx, :);
    matrix(min_idx, :) = [];
    cur_line = cur_line + 1;
end
    


