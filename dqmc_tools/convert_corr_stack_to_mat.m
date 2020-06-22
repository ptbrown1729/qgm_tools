function corr_mat = convert_corr_stack_to_mat(corr_struct_stack)
    inds = [];
    corr_regex = 'point__dx__(\d+)__dy__(\d+)';
    struct_fields = fields(corr_struct_stack);
    
    for ii = 1:length(struct_fields)
        m = regexp(struct_fields{ii}, corr_regex, 'tokens');
        inds = vertcat(inds, [str2num(m{1}{1}), str2num(m{1}{2})]);
    end

    mat_size = max(inds(:)) + 1;
    corr_mat = zeros(length(corr_struct_stack), mat_size, mat_size, 2);
    for ii = 1:length(inds)
        xx = inds(ii, 1) + 1; % matlab indexing starts at 1 instead of 0
        yy = inds(ii, 2) + 1;
        corr_mat(:, yy, xx, :) = vertcat(corr_struct_stack.(struct_fields{ii}));
        corr_mat(:, xx, yy, :) = corr_mat(:, yy, xx, :);
    end
end