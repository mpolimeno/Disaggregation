function  [cubes_that_break,where] = GetCubesWhoseFaceDetaches(face_detached,...
                                                               internal_faces_and_cubes_index_array,...
                                                               internal_faces_and_cubes_index_array_no_double_counting)

    % Index of the face that detaches
    face_index = find(internal_faces_and_cubes_index_array(:,1)==face_detached);

    % This gives the cubes that break
    cubes_that_break = internal_faces_and_cubes_index_array(face_index,2:3);

    % Sort third column based on second column to feed columns to Graph (which sorts automatically)
    % This avoids indexing issues
    internal_faces_and_cubes_index_array_no_double_counting = sortrows(internal_faces_and_cubes_index_array_no_double_counting,[2 3])
    % get index of which cubes break in the unrepeated array
    for ii=1:length(internal_faces_and_cubes_index_array_no_double_counting)
        internal_faces_and_cubes_index_array_no_double_counting(ii,2:3)
        cubes_that_break(1:2)
        cube_bool = (internal_faces_and_cubes_index_array_no_double_counting(ii,2:3)==cubes_that_break(1:2));
        if sum(cube_bool)==2
            where = ii;
            break;
        end
    end

end
