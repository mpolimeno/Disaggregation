function force_acting_on_cube = ComputeForceActingOnCube(xc,drag,torque,flow)
    
    % The scaling argument about the velocity and the force/torque does not seem to be what we want
    % Here is the proposed new implementation
    % If we have only a drag acting on the object, then the force on each cube will simply be drag/NC
    % If we have only a torque acting on the object, then the force on each cube will be given by torque\cross(x_i-x_cm)/||(x_i-x_cm)||^2

    NC = size(xc,1);
    if size(xc,1)==1
        cm = xc; % for one cube only
    else
        cm = mean(xc); % for multiple cubes
    end
    xc_shift = xc - cm;

    % In the case of extensional flow the aggregate will not move, so the velocity of each cube will be identically zero
    % If I am returning a zero force on each cube, this first condition is actually not necessary.
    % You can even say that even calling the function in the first place is not necessary, but should probably be done for logic consistency
    force_acting_on_cube = zeros(NC,3);
    if flow==3 || flow==4
        force_acting_on_cube = zeros(NC,3);
    elseif flow==1
        for cube_index=1:NC
            force_acting_on_cube(cube_index,:) = drag'/NC;
        end
    else
        for cube_index=1:NC
            x_vec = (xc_shift(cube_index,:));
            force_acting_on_cube(cube_index,:) = cross(torque',x_vec)/(norm(x_vec))^2;
        end
    end
end
