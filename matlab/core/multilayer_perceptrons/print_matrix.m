function print_matrix(fileID, matrix)
    
    format_str = ['  %', num2str(ceil(log10(max(max(matrix))))), '.0f'];
    for r = 1:size(matrix, 1)
        for c = 1:size(matrix, 2)
            fprintf(fileID, format_str, matrix(r,c));
        end
        fprintf(fileID, '\n');
    end
end