function d2_indices = get_d2_indices(initial_index)
    %get all symmetric indices assuming D4 (i.e. square)
    %symmetry.
    %
    % initial_index: represents a coordinate
    
    d2_indices = zeros(4, 2);
    d2_indices(1,:) = initial_index;
    d2_indices(2,:) = -initial_index;
    d2_indices(3,:) = [-initial_index(1), initial_index(2)];
    d2_indices(4,:) = [initial_index(1), -initial_index(2)];
    
    d2_indices = unique(d2_indices, 'rows');
    
end