function PlotInternalStresses(center_of_cubes,...
                              internal_and_external_stresses,...
                              total_number_of_cubes,...
                              figure_number,...
                              center_of_external_faces,...
                              total_number_of_external_faces,...
                              cubes_that_break)
                        
    
    %{
     This function plots the internal stresses for a fractal aggregate made
     of a given number of cubes.
     Inputs:
        center_of_cubes {array}: the coordinates of the center (x,y,z) of
        each cube in the given aggregate;
        internal_and_external_stresses {array}: values of the stresses computed on all
        the faces;
    %}

    figure(figure_number)
    clf
    xlabel("x")
    ylabel("y")
    zlabel("z")
    
    number_of_faces_in_a_single_cube = 6;
    norm_of_all_forces = zeros(number_of_faces_in_a_single_cube*total_number_of_cubes,1);
    for i=1:(number_of_faces_in_a_single_cube*total_number_of_cubes)
        norm_of_all_forces(i) = norm(internal_and_external_stresses(:,i));
    end

    max_stress = max(norm_of_all_forces);
    min_stress = min(norm_of_all_forces);

    range_of_stresses = max_stress-min_stress;

    max_value_of_stress = 0;
    for i=1:total_number_of_cubes

        % face 1
        center_of_face        = [center_of_cubes(i,1)+1 center_of_cubes(i,2) center_of_cubes(i,3)]; % face center
        x_vertices_of_face_s  = [center_of_cubes(i,1)+1 center_of_cubes(i,1)+1 center_of_cubes(i,1)+1 center_of_cubes(i,1)+1];
        y_vertices_of_face_s  = [center_of_cubes(i,2)+1 center_of_cubes(i,2)+1 center_of_cubes(i,2)-1 center_of_cubes(i,2)-1];
        z_vertices_of_face_s  = [center_of_cubes(i,3)+1 center_of_cubes(i,3)-1 center_of_cubes(i,3)-1 center_of_cubes(i,3)+1];
        stress_on_face = (norm(internal_and_external_stresses(:,6*(i-1)+1))-min_stress)/range_of_stresses;
        if (stress_on_face > max_value_of_stress)
            max_value_of_stress = stress_on_face;
        end
        this_is_an_internal_face = true; % assume we are inside
        for j=1:total_number_of_external_faces
            if (center_of_face==center_of_external_faces(j,:)) % we are on the outside
                this_is_an_internal_face = false;
            end 
            if (this_is_an_internal_face==false)
                patch(x_vertices_of_face_s,...
                      y_vertices_of_face_s,...
                      z_vertices_of_face_s,...
                      [1 1 1],...
                      'FaceAlpha',0)
            end
        end
        if (this_is_an_internal_face==true)
            patch(x_vertices_of_face_s,...
                  y_vertices_of_face_s,...
                  z_vertices_of_face_s,...
                  [1-stress_on_face 1-stress_on_face 1])
            %cube_where_max_stress_is_located = i;
            face_where_max_stress_is_located = 1;
        end

        % face 2
        center_of_face        = [center_of_cubes(i,1) center_of_cubes(i,2)+1 center_of_cubes(i,3)]; % face center
        x_vertices_of_face_s  = [center_of_cubes(i,1)+1 center_of_cubes(i,1)+1 center_of_cubes(i,1)-1 center_of_cubes(i,1)-1];
        y_vertices_of_face_s  = [center_of_cubes(i,2)+1 center_of_cubes(i,2)+1 center_of_cubes(i,2)+1 center_of_cubes(i,2)+1];
        z_vertices_of_face_s  = [center_of_cubes(i,3)+1 center_of_cubes(i,3)-1 center_of_cubes(i,3)-1 center_of_cubes(i,3)+1];
        stress_on_face = (norm(internal_and_external_stresses(:,6*(i-1)+2))-min_stress)/range_of_stresses;
        if (stress_on_face > max_value_of_stress)
            max_value_of_stress = stress_on_face;
        end
        this_is_an_internal_face = true; % assume we are inside
        for j=1:total_number_of_external_faces
            if (center_of_face==center_of_external_faces(j,:)) % we are on the outside
                this_is_an_internal_face = false;
            end  
            if (this_is_an_internal_face==false)
                patch(x_vertices_of_face_s,...
                      y_vertices_of_face_s,...
                      z_vertices_of_face_s,...
                      [1 1 1],...
                      'FaceAlpha',0)
            end
        end
        if (this_is_an_internal_face==true)
            patch(x_vertices_of_face_s,...
                  y_vertices_of_face_s,...
                  z_vertices_of_face_s,...
                  [1-stress_on_face 1-stress_on_face 1])
            %cube_where_max_stress_is_located = i;
            face_where_max_stress_is_located = 2;
        end

        % face 3
        center_of_face        = [center_of_cubes(i,1) center_of_cubes(i,2) center_of_cubes(i,3)+1]; % face center
        x_vertices_of_face_s  = [center_of_cubes(i,1)+1 center_of_cubes(i,1)-1 center_of_cubes(i,1)-1 center_of_cubes(i,1)+1];
        y_vertices_of_face_s  = [center_of_cubes(i,2)+1 center_of_cubes(i,2)+1 center_of_cubes(i,2)-1 center_of_cubes(i,2)-1];
        z_vertices_of_face_s  = [center_of_cubes(i,3)+1 center_of_cubes(i,3)+1 center_of_cubes(i,3)+1 center_of_cubes(i,3)+1];
        stress_on_face = (norm(internal_and_external_stresses(:,6*(i-1)+3))-min_stress)/range_of_stresses;
        if (stress_on_face > max_value_of_stress)
            max_value_of_stress = stress_on_face;
        end
        this_is_an_internal_face = true; % assume we are inside
        for j=1:total_number_of_external_faces
            if (center_of_face==center_of_external_faces(j,:)) % we are on the outside
                this_is_an_internal_face = false;
            end
            if (this_is_an_internal_face==false)
                patch(x_vertices_of_face_s,...
                      y_vertices_of_face_s,...
                      z_vertices_of_face_s,...
                      [1 1 1],...
                      'FaceAlpha',0)
            end
        end
        if this_is_an_internal_face==true
            patch(x_vertices_of_face_s,...
                  y_vertices_of_face_s,...
                  z_vertices_of_face_s,...
                  [1-stress_on_face 1-stress_on_face 1])
            %cube_where_max_stress_is_located = i;
            face_where_max_stress_is_located = 3;
        end

        % face 4
        center_of_face        = [center_of_cubes(i,1) center_of_cubes(i,2) center_of_cubes(i,3)-1]; % face center
        x_vertices_of_face_s  = [center_of_cubes(i,1)+1 center_of_cubes(i,1)-1 center_of_cubes(i,1)-1 center_of_cubes(i,1)+1];
        y_vertices_of_face_s  = [center_of_cubes(i,2)+1 center_of_cubes(i,2)+1 center_of_cubes(i,2)-1 center_of_cubes(i,2)-1];
        z_vertices_of_face_s  = [center_of_cubes(i,3)-1 center_of_cubes(i,3)-1 center_of_cubes(i,3)-1 center_of_cubes(i,3)-1];
        stress_on_face = (norm(internal_and_external_stresses(:,6*(i-1)+4))-min_stress)/range_of_stresses;
        if (stress_on_face > max_value_of_stress)
            max_value_of_stress = stress_on_face;
        end
        this_is_an_internal_face = true; % assume we are inside
        for j=1:total_number_of_external_faces
            if (center_of_face==center_of_external_faces(j,:)) % we are on the outside
                this_is_an_internal_face = false;
            end
            if (this_is_an_internal_face==false)
                patch(x_vertices_of_face_s,...
                      y_vertices_of_face_s,...
                      z_vertices_of_face_s,...
                      [1 1 1],...
                      'FaceAlpha',0)
            end
        end
        if this_is_an_internal_face==true
            patch(x_vertices_of_face_s,...
                  y_vertices_of_face_s,...
                  z_vertices_of_face_s,...
                  [1-stress_on_face 1-stress_on_face 1])
            %cube_where_max_stress_is_located = i;
            face_where_max_stress_is_located = 4;
        end

        % face 5
        center_of_face        = [center_of_cubes(i,1) center_of_cubes(i,2)-1 center_of_cubes(i,3)]; % face center
        x_vertices_of_face_s  = [center_of_cubes(i,1)+1 center_of_cubes(i,1)+1 center_of_cubes(i,1)-1 center_of_cubes(i,1)-1];
        y_vertices_of_face_s  = [center_of_cubes(i,2)-1 center_of_cubes(i,2)-1 center_of_cubes(i,2)-1 center_of_cubes(i,2)-1];
        z_vertices_of_face_s  = [center_of_cubes(i,3)+1 center_of_cubes(i,3)-1 center_of_cubes(i,3)-1 center_of_cubes(i,3)+1];
        stress_on_face = (norm(internal_and_external_stresses(:,6*(i-1)+5))-min_stress)/range_of_stresses;
        if (stress_on_face > max_value_of_stress)
            max_value_of_stress = stress_on_face;
        end
        this_is_an_internal_face = true; % assume we are inside
        for j=1:total_number_of_external_faces
            if (center_of_face==center_of_external_faces(j,:)) % we are on the outside
                this_is_an_internal_face = false;
            end
            if (this_is_an_internal_face==false)
                patch(x_vertices_of_face_s,...
                      y_vertices_of_face_s,...
                      z_vertices_of_face_s,...
                      [1 1 1],...
                      'FaceAlpha',0)
            end
        end
        if this_is_an_internal_face==true
            patch(x_vertices_of_face_s,...
                  y_vertices_of_face_s,...
                  z_vertices_of_face_s,...
                  [1 1 1])
            %cube_where_max_stress_is_located = i;
            face_where_max_stress_is_located = 5;
        end

        % face 6
        center_of_face        = [center_of_cubes(i,1)-1 center_of_cubes(i,2) center_of_cubes(i,3)]; % face center
        x_vertices_of_face_s  = [center_of_cubes(i,1)-1 center_of_cubes(i,1)-1 center_of_cubes(i,1)-1 center_of_cubes(i,1)-1];
        y_vertices_of_face_s  = [center_of_cubes(i,2)+1 center_of_cubes(i,2)+1 center_of_cubes(i,2)-1 center_of_cubes(i,2)-1];
        z_vertices_of_face_s  = [center_of_cubes(i,3)+1 center_of_cubes(i,3)-1 center_of_cubes(i,3)-1 center_of_cubes(i,3)+1];
        stress_on_face = (norm(internal_and_external_stresses(:,6*(i-1)+6))-min_stress)/range_of_stresses;
        if (stress_on_face > max_value_of_stress)
            max_value_of_stress = stress_on_face;
        end
        this_is_an_internal_face = true; % assume we are inside
        for j=1:total_number_of_external_faces
            if (center_of_face==center_of_external_faces(j,:)) % we are on the outside
                this_is_an_internal_face = false;
            end
            if (this_is_an_internal_face==false)
                patch(x_vertices_of_face_s,...
                      y_vertices_of_face_s,...
                      z_vertices_of_face_s,...
                      [1 1 1],...
                      'FaceAlpha',0)
            end
        end
        if this_is_an_internal_face==true
            patch(x_vertices_of_face_s,...
                  y_vertices_of_face_s,...
                  z_vertices_of_face_s,...
                  [1-stress_on_face 1-stress_on_face 1])
            %cube_where_max_stress_is_located = i;
            face_where_max_stress_is_located = 6;
        end

    end

    % Max stress will be on a red face
    if (face_where_max_stress_is_located==1)
        cube_where_max_stress_is_located = cubes_that_break(1);
        x_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,1)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)+1];
                            
        y_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,2)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)-1];
                            
        z_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,3)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)+1];
                            
        patch(x_vertices_of_face_s,y_vertices_of_face_s,z_vertices_of_face_s,[1 0 0]) % this is red
    end
    
    if (face_where_max_stress_is_located==2)
        cube_where_max_stress_is_located = cubes_that_break(1);
        x_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,1)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)-1];
                            
        y_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,2)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)+1];
                            
        z_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,3)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)+1];
                            
        patch(x_vertices_of_face_s,y_vertices_of_face_s,z_vertices_of_face_s,[1 0 0]) % this is red
    end
    
    if (face_where_max_stress_is_located==3)
        cube_where_max_stress_is_located = cubes_that_break(1);
        x_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,1)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)+1];
                            
        y_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,2)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)-1];
                            
        z_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,3)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)+1];
                            
        patch(x_vertices_of_face_s,y_vertices_of_face_s,z_vertices_of_face_s,[1 0 0]) % this is red
    end
    
    if (face_where_max_stress_is_located==4)
        cube_where_max_stress_is_located = cubes_that_break(1);
        x_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,1)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)+1];
                            
        y_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,2)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)-1];
                            
        z_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,3)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)-1];
                            
        patch(x_vertices_of_face_s,y_vertices_of_face_s,z_vertices_of_face_s,[1 0 0]) % this is red
    end
    
    if (face_where_max_stress_is_located==5)
        cube_where_max_stress_is_located = cubes_that_break(1);
        x_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,1)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)-1];
                            
        y_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,2)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)-1];
                            
        z_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,3)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)+1];
                            
        patch(x_vertices_of_face_s,y_vertices_of_face_s,z_vertices_of_face_s,[1 0 0]) % this is red
    end
    
    if (face_where_max_stress_is_located==6)
        cube_where_max_stress_is_located = cubes_that_break(1);
        x_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,1)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,1)-1];
                            
        y_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,2)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,2)-1];
                            
        z_vertices_of_face_s = [center_of_cubes(cube_where_max_stress_is_located,3)+1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)-1,...
                                center_of_cubes(cube_where_max_stress_is_located,3)+1];
                            
        patch(x_vertices_of_face_s,y_vertices_of_face_s,z_vertices_of_face_s,[1 0 0]) % this is red
    end

    axis equal
end

