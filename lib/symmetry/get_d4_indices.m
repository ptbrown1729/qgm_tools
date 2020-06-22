function d4_indices = get_d4_indices(initial_index)
    %get all symmetric indices assuming D4 (i.e. square)
    %symmetry.
    %
    % initial_index: represents a coordinate
    
    d4_indices = zeros(8, 2);
    d4_indices(1,:) = initial_index;
    d4_indices(2,:) = -initial_index;
    d4_indices(3,:) = [initial_index(2), initial_index(1)];
    d4_indices(4,:) = -d4_indices(3,:);
    d4_indices(5,:) = [-initial_index(1), initial_index(2)];
    d4_indices(6,:) = [initial_index(1), -initial_index(2)];
    d4_indices(7,:) = [-initial_index(2), initial_index(1)];
    d4_indices(8,:) = [initial_index(2), -initial_index(1)];
    
    d4_indices = unique(d4_indices, 'rows');
    
end