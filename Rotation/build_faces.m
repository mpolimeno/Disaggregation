function [center_of_external_faces, normal_direction_of_external_faces,...
          orientation_of_external_faces,...
          number_of_external_faces] = build_faces(xc, total_number_of_cubes)
    % builds the surface of the object make of cubes of side 2 with centers at location xc and number NC
    % This files builds the face centers, normal directions, and orientations of a volume
    % made of cubes sharing at least one face.
    % Thankfully, the order does not matter.

    % clear('posint')

    % NC = 6; % number of cubes
    % xc(i,:) % the position of the centers, of the form (2*j,2*k,2*m), for j, k, m elements of Z
    % xc(1,:) = [0 0 0];
    % xc(2,:) = [2 0 0];
    % xc(3,:) = [4 0 0];
    % xc(4,:) = [2 2 0];
    % xc(5,:) = [0 2 0];
    % xc(6,:) = [0 0 2];


    % base cube with one point per face
    position_of_faces_of_base_cube(1,:) = [1,0,0];
    position_of_faces_of_base_cube(2,:) = [0,1,0];
    position_of_faces_of_base_cube(3,:) = [0,0,1];
    position_of_faces_of_base_cube(4,:) = [0,0,-1];
    position_of_faces_of_base_cube(5,:) = [0,-1,0];
    position_of_faces_of_base_cube(6,:) = [-1,0,0];
    
    % Normal directions and orientation of each face of the base cube
    normal_direction_of_each_face(1)    = 1;
    orientation_of_each_normal(1)       = 1;
    
    normal_direction_of_each_face(2)    = 2;
    orientation_of_each_normal(2)       = 1;
    
    normal_direction_of_each_face(3)    = 3;
    orientation_of_each_normal(3)       = 1;
    
    normal_direction_of_each_face(4)    = 3;
    orientation_of_each_normal(4)       = -1;
    
    normal_direction_of_each_face(5)    = 2;
    orientation_of_each_normal(5)       = -1;
    
    normal_direction_of_each_face(6)    = 1;
    orientation_of_each_normal(6)       = -1;


    centers_of_external_faces_list          = []; % to hold face centers of all cubes
    normal_direction_of_external_faces_list = []; % to hold normal directions of all cubes
    orientation_of_external_faces_list      = []; % to hold normal orirentations of all cubes
    
    index_of_faces_to_be_removed_from_cube  = [];


    number_of_faces_in_one_cube = 6;

    for ii=1:total_number_of_cubes
      % add all 6 face centers
        for kk = 1:number_of_faces_in_one_cube
            centers_of_external_faces_list = [centers_of_external_faces_list; xc(ii,:)+position_of_faces_of_base_cube(kk,:)];
            normal_direction_of_external_faces_list = [normal_direction_of_external_faces_list; normal_direction_of_each_face(kk)];
            orientation_of_external_faces_list = [orientation_of_external_faces_list; orientation_of_each_normal(kk)];
        end
        % check if there is an cube adjacent cube already
        for jj=1:ii-1
            for kk = 1:number_of_faces_in_one_cube
                if (xc(ii,:) + position_of_faces_of_base_cube(kk,:)*2 == xc(jj,:))
                    index_of_faces_to_be_removed_from_cube = [index_of_faces_to_be_removed_from_cube; ii,kk]; % remove that face from the new cube
                    index_of_faces_to_be_removed_from_cube = [index_of_faces_to_be_removed_from_cube; jj,7-kk]; % remove corresponding face from the old cube
                end
            end
        end
    end

    % Now remove all the useless faces
    number_of_faces_to_remove = length(index_of_faces_to_be_removed_from_cube);
    index_of_faces_to_keep = [];
    for jj=1:(number_of_faces_in_one_cube*total_number_of_cubes)
        remove_this_face = false;
        for ii=1:number_of_faces_to_remove
            current_index = number_of_faces_in_one_cube*(index_of_faces_to_be_removed_from_cube(ii,1)-1)...
                   +index_of_faces_to_be_removed_from_cube(ii,2);
            if (jj==current_index)
                remove_this_face = true;
            end
        end
        if (remove_this_face==false)
            index_of_faces_to_keep = [index_of_faces_to_keep jj];
        end
    end

    center_of_external_faces = centers_of_external_faces_list(index_of_faces_to_keep,:);
    normal_direction_of_external_faces = normal_direction_of_external_faces_list(index_of_faces_to_keep);
    orientation_of_external_faces = orientation_of_external_faces_list(index_of_faces_to_keep);

    number_of_external_faces = length(center_of_external_faces);

end







