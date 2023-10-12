function PlotInternalStresses(figure_number,...
                              center_of_external_faces,...
                              finalndir,...
                              stress_and_faces,...
                              input_array)
                        
    
    %{
     This function plots the internal stresses for a fractal aggregate made
     of a given number of cubes.
     Inputs:
        figure_number {int}: User-input figure number
        center_of_external_faces {array}: the coordinates of the center (x,y,z) of
        each external face in the given aggregate;
        finalndir {array}: direction of normal vector in each external face
        input_array {array}: contains indices of internal faces, their cubes, 
                             their associated direction and their associated stress 
    %}

    figure(figure_number)
    clf
    xlabel('x',"Interpreter","LaTex")
    ylabel('y',"Interpreter","LaTex")
    zlabel('z',"Interpreter","LaTex")
    set(gca,"FontSize",25)
    
    % add path to colormap
    addpath(genpath('../matplotlib'));
    c = colormap(flipud(viridis(size(stress_and_faces,1))));
    my_colormap = ([c stress_and_faces(:,2)]);
    colorbar

    center_of_internal_faces = input_array(:,1:3);
    finalndir_internal = input_array(:,4);
    
    total_number_of_external_faces = size(center_of_external_faces,1);
    total_number_of_internal_faces = size(center_of_internal_faces,1);
    % let us take care of the external faces first
    % they will all be transparent
    for external=1:total_number_of_external_faces
        if (finalndir(external) == 1)
            x_vertices_of_face_s = [center_of_external_faces(external,1) center_of_external_faces(external,1) center_of_external_faces(external,1) center_of_external_faces(external,1) ];
            y_vertices_of_face_s = [center_of_external_faces(external,2)+1 center_of_external_faces(external,2)+1 center_of_external_faces(external,2)-1 center_of_external_faces(external,2)-1 ];
            z_vertices_of_face_s = [center_of_external_faces(external,3)+1 center_of_external_faces(external,3)-1 center_of_external_faces(external,3)-1 center_of_external_faces(external,3)+1 ];
            hold on
        end
        if (finalndir(external) == 2)
            x_vertices_of_face_s = [center_of_external_faces(external,1)+1 center_of_external_faces(external,1)+1 center_of_external_faces(external,1)-1 center_of_external_faces(external,1)-1 ];
            y_vertices_of_face_s = [center_of_external_faces(external,2) center_of_external_faces(external,2) center_of_external_faces(external,2) center_of_external_faces(external,2) ];
            z_vertices_of_face_s = [center_of_external_faces(external,3)+1 center_of_external_faces(external,3)-1 center_of_external_faces(external,3)-1 center_of_external_faces(external,3)+1 ];
            hold on
        end
        if (finalndir(external) == 3)
            x_vertices_of_face_s = [center_of_external_faces(external,1)+1 center_of_external_faces(external,1)-1 center_of_external_faces(external,1)-1 center_of_external_faces(external,1)+1 ];
            y_vertices_of_face_s = [center_of_external_faces(external,2)+1 center_of_external_faces(external,2)+1 center_of_external_faces(external,2)-1 center_of_external_faces(external,2)-1 ];
            z_vertices_of_face_s = [center_of_external_faces(external,3) center_of_external_faces(external,3) center_of_external_faces(external,3) center_of_external_faces(external,3) ];
            hold on
        end
        patch(x_vertices_of_face_s,...
              y_vertices_of_face_s,...
              z_vertices_of_face_s,...
              [1 1 1],...
              'FaceAlpha',0)
    end
    
    hold on

    for internal=1:total_number_of_internal_faces
        if (finalndir_internal(internal) == 1)
            x_vertices_of_face_s = [center_of_internal_faces(internal,1) center_of_internal_faces(internal,1) center_of_internal_faces(internal,1) center_of_internal_faces(internal,1) ];
            y_vertices_of_face_s = [center_of_internal_faces(internal,2)+1 center_of_internal_faces(internal,2)+1 center_of_internal_faces(internal,2)-1 center_of_internal_faces(internal,2)-1 ];
            z_vertices_of_face_s = [center_of_internal_faces(internal,3)+1 center_of_internal_faces(internal,3)-1 center_of_internal_faces(internal,3)-1 center_of_internal_faces(internal,3)+1 ];
            hold on
        end
        if (finalndir_internal(internal) == 2)
            x_vertices_of_face_s = [center_of_internal_faces(internal,1)+1 center_of_internal_faces(internal,1)+1 center_of_internal_faces(internal,1)-1 center_of_internal_faces(internal,1)-1 ];
            y_vertices_of_face_s = [center_of_internal_faces(internal,2) center_of_internal_faces(internal,2) center_of_internal_faces(internal,2) center_of_internal_faces(internal,2) ];
            z_vertices_of_face_s = [center_of_internal_faces(internal,3)+1 center_of_internal_faces(internal,3)-1 center_of_internal_faces(internal,3)-1 center_of_internal_faces(internal,3)+1 ];
            hold on
        end
        if (finalndir_internal(internal) == 3)
            x_vertices_of_face_s = [center_of_internal_faces(internal,1)+1 center_of_internal_faces(internal,1)-1 center_of_internal_faces(internal,1)-1 center_of_internal_faces(internal,1)+1 ];
            y_vertices_of_face_s = [center_of_internal_faces(internal,2)+1 center_of_internal_faces(internal,2)+1 center_of_internal_faces(internal,2)-1 center_of_internal_faces(internal,2)-1 ];
            z_vertices_of_face_s = [center_of_internal_faces(internal,3) center_of_internal_faces(internal,3) center_of_internal_faces(internal,3) center_of_internal_faces(internal,3) ];
            hold on
        end
        face_index = stress_and_faces(internal,2);
        what_color = find(my_colormap(:,4)==face_index);
        patch(x_vertices_of_face_s,...
              y_vertices_of_face_s,...
              z_vertices_of_face_s,...
              [my_colormap(what_color,1)...
               my_colormap(what_color,2)...
               my_colormap(what_color,3)])
    end
    axis equal
    az = 121;
    el = 31;
    view(az,el)
end