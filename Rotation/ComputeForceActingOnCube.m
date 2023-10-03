function force_acting_on_cube = ComputeForceActingOnCube(xc,U_vec,Omega_vec,drag,torque,flow)
    NC = size(xc,1);
    unit_drag   = drag/NC;
    unit_torque = torque/NC;
    if size(xc,1)==1
        cm = xc; % for one cube only
    else
        cm = mean(xc); % for multiple cubes
    end

    % I think I can just shift everything back just for this
    xc_shift = xc - cm;

    % I need to make this into a function
    % that function will have the shifted coordinate system
    % so that I can compute the velocity and the force on the cubes
    % then I shift back to integers to get the correct indices of the faces
    velocity_of_cube     = [];
    force_acting_on_cube = zeros(NC,3);
    x_vec_list = [];

    % In the case of extensional flow the aggregate will not move, so the velocity of each cube will be identically zero
    % If I am returning a zero force on each cube, this first condition is actually not necessary.
    % You can even say that even calling the function in the first place is not necessary, but should probably be done for logic consistency
    if flow==3 || flow==4
        u_i = zeros(NC,3);
        velocity_of_cube = [velocity_of_cube;u_i];
    else
        for jj=1:NC
            x_vec = (xc_shift(jj,:));
            if flow==1
                u_i = U_vec';
            elseif flow==2
                u_i = cross(Omega_vec',x_vec);
            end
            velocity_of_cube = [velocity_of_cube;u_i];
            x_vec_list = [x_vec_list;x_vec];
        end
    end

    if flow==1 || flow==2
        proportionality_constant_for_translation = zeros(NC,1);
        proportionality_constant_for_rotation    = zeros(NC,1);
        if round(norm(drag))~=0 % We are in the case in which we impose a drag, but not a torque. Note: we always impose either or
            for cube_index = 1:NC
                proportionality_constant_for_translation(cube_index) = norm(velocity_of_cube(cube_index,:))/norm(unit_drag);
                if norm(proportionality_constant_for_translation(cube_index))==0
                    force_acting_on_cube(cube_index,1:3) = [0,0,0];
                else
                    force_acting_on_cube(cube_index,1:3) = -velocity_of_cube(cube_index,:)/proportionality_constant_for_translation(cube_index);
                end
            end
        else
            for cube_index=1:NC
                force_acting_on_cube(cube_index,:) = cross(torque',x_vec_list(cube_index,:))/(norm(x_vec_list(cube_index,:)))^2;
%                 proportionality_constant_for_rotation(cube_index) = norm(velocity_of_cube(cube_index,:))/norm(unit_torque);
%                 if norm(proportionality_constant_for_rotation(cube_index))==0
%                     force_acting_on_cube(cube_index,1:3) = [0,0,0];
%                 else
%                     force_acting_on_cube(cube_index,1:3) = -velocity_of_cube(cube_index,:)/proportionality_constant_for_rotation(cube_index);
%                 end
            end
        end
    end
end
