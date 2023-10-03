function internal_faces = FindInternalFaces(xc,NC)

    % This function returns the index of the internal faces of an aggregate
    % with NC cubes, whose center is at xc
    % if size(xc,1)==1
    %     cm = xc; % for one cube only
    % else
    %     cm = mean(xc); % for multiple cubes
    % end

    base_cube(1,:) = [1,0,0];
    base_cube(2,:) = [0,1,0];
    base_cube(3,:) = [0,0,1];
    base_cube(4,:) = [0,0,-1];
    base_cube(5,:) = [0,-1,0];
    base_cube(6,:) = [-1,0,0];

    % index_of_any_face is an array of flags:
    % If it is 2, it is an outer face,
    % if it is 1, it is an inner face already dealt with
    % if it is 0, it is a new inner face
    index_of_unknown_face = 0;
    index_for_outer_faces = 0; % index for outer faces

    % flags for each face
    index_of_any_face = zeros(1,NC*6); % to check if we have an outer face or an inner face

    position_of_neighbor_list = [];
    internal_faces_and_cubes_index_array = [];
    internal_faces_and_cubes_index_array_no_double_counting = [];
    for i=1:NC
        for k = 1:6
            index_of_unknown_face = index_of_unknown_face + 3; % update index of faces, indexes of unknowns
            position_of_neighbor = xc(i,:)+2*base_cube(k,:);
            if (index_of_any_face(index_of_unknown_face/3) == 0) % never been to that face
                for j=1:NC % check if it is an outer face
                    if ( position_of_neighbor == xc(j,:)) % We have a neighbor
                        position_of_neighbor_list = [position_of_neighbor_list;position_of_neighbor];
                        index_of_any_face(index_of_unknown_face/3) = 1; % this face was visited and is an inner face
                        internal_faces_and_cubes_index_array_no_double_counting = [internal_faces_and_cubes_index_array_no_double_counting; (index_of_unknown_face/3) i j];
                        internal_faces_and_cubes_index_array = [internal_faces_and_cubes_index_array; (index_of_unknown_face/3) i j];
                        index_of_unknown_face2 = (6*(j-1)+7-k)*3; % index of the corresponding inner face we just found
                        index_of_any_face(index_of_unknown_face2/3) = 1; % the corresponding face was visited too
                        internal_faces_and_cubes_index_array = [internal_faces_and_cubes_index_array; (index_of_unknown_face2/3) i j];
                        break; % because we found our neighbor, we can stop looking not sure if that works ??
                   end
               end
              if (index_of_any_face(index_of_unknown_face/3)==0) % We have no neighbor
                  %index_of_unknown_face
                index_of_any_face(index_of_unknown_face/3) = 2; % this face is an outer face
                index_for_outer_faces = index_for_outer_faces + 1;
              end
            end
        end
    end

    % this will give me the index of the internal faces
    internal_faces = find(index_of_any_face==1);
    % index_of_any_face is an array of flags -> if it stores 1, that is an internal face

end
