function PlotMaximumInternalStress(figure_number,...
                              center_of_external_faces,...
                              finalndir,...
                              stress_and_faces,...
                              input_array,...
                              face_detached)
                        
    
    %{
     This function highlighs the location of the maximum internal stress
     for a given fractal aggregate
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

    center_of_internal_faces = input_array(:,1:3);
    finalndir_internal = input_array(:,4);
    
    total_number_of_external_faces = size(center_of_external_faces,1);
    total_number_of_internal_faces = size(center_of_internal_faces,1);
    % let us take care of the external faces first
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
              'FaceAlpha',0,...
              "LineWidth",1)
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
        if face_index==face_detached
            what_color = [1 0 0]; % red
            patch(x_vertices_of_face_s,...
                  y_vertices_of_face_s,...
                  z_vertices_of_face_s,...
                  [what_color(1)...
                   what_color(2)...
                   what_color(3)],...
                   "LineWidth",2)
        else
            what_color = [1 1 1];
            patch(x_vertices_of_face_s,...
                  y_vertices_of_face_s,...
                  z_vertices_of_face_s,...
                  [what_color(1)...
                   what_color(2)...
                   what_color(3)],...
                  "FaceAlpha",0,...
                  "LineWidth",1) % transparent
        end
    end
    axis equal
    az = 121;
    el = 31;
    view(az,el)
end