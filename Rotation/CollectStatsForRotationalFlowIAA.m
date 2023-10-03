NC_vec = [10,25,50,100,150,200,250,300,400];% number of cubes
SEED_vec = 1:50;

count = 0;
loop_count = 0;
for nc=1:length(NC_vec)
    count = count + 1;
    for seed=1:length(SEED_vec)
        NC = NC_vec(nc);
        SEED = SEED_vec(seed);
        %xc = DLA_3D(NC,SEED);
        
        position_data = load(sprintf("~/Desktop/Disaggregation/IAA/Stats/SEED_%i_ORIGINAL_Agg_NC_%i_BROKEN_FLOW_1.mat",SEED,NC));
        xc = position_data.xc;

        if size(xc,1)==1
            cm = xc; % for one cube only
        else
            cm = mean(xc); % for multiple cubes
        end

        % Now compute where are the faces and what are their normals and orientations
        [finalposint, finalndir, finalori,Nf] = build_faces(xc, NC);
        % Nf is the number of faces
        % finalposint is the location of the center of each face
        % finalndir is the direction of the normal
        % finalori is the orientation of each face
        % Nf is the number of faces

    %%%%%%%%%%%%%%%% SELECT FLOW %%%%%%%%%%%%%%%%%%%

    % all 4 of these flows recover the correct velocities
    % flow = 1 -> translation
    % flow = 2 -> rotation
    % flow = 3 -> extension 
    % flow = 4 -> shear
    flow = 2;

    if flow~=1 && flow~=2 && flow~=3 && flow~=4
        msg = "Flow must be either 1, 2, 3, 4";
        error(msg);
    end
    % define matrix for extensional flow
    % from here it will also get passed to Double Layer solver
    % if the flow is of other type, pass empty vector
    if flow==3
       M   = [-1,0,0;0,1,0;0,0,0];
    elseif flow==4
        M = [0,1,0;1,0,0;0,0,0]; % for shear
    else
        M = [];
    end

    % based on different flows, select appropriate codes to get drag and torque
    % to be fed to the new code with double layer
    if flow==1
        msg = "Translational flow selected";
        disp(msg);
        U_translation = [0;0;1]; % translational velocity
        [forceout,drag,torque] = fractal_bi_stokes_force(xc,finalposint,finalndir,finalori,U_translation,Nf);
        % fixing drag, torque and background flow
        drag_in         = drag;
        torque_in       = [0;0;0];
        background_flow = [0;0;0];
    elseif flow==2
        msg = "Rotational flow selected";
        disp(msg);
        Omega_rotation = [0;0;1]; % rotational velocity
        [forceout,drag,torque] = fractal_bi_stokes_force_rot(xc,finalposint,finalndir,finalori,Omega_rotation,Nf);
        % fixing drag, torque and background flow
        drag_in         = [0;0;0];
        torque_in       = torque;
        background_flow = [0;0;0];
    elseif flow==3
        msg = "Extensional flow selected";
        disp(msg);
        % fixing drag, torque and background flow
        drag_in            = [0;0;0];
        torque_in          = [0;0;0];
        background_flow    = [0;0;0];
    elseif flow==4
        msg = "Shear flow selected";
        disp(msg);
        % fixing drag, torque and background flow
        drag_in         = [0;0;0];
        torque_in       = [0;0;0];
        background_flow = 0.5*[0;0;-1];
    else
        error_message = "ERROR: FLOW NOT READY YET";
        error(error_message);
    end

    % Shift origin of system to center of mass of aggregate
    xc_shift = xc - cm;
    finalposint_shift = finalposint - cm;

    % Compute External Stresses and Velocities for solid body motion
    [LHS,sol,stress_outer,U_vec,Omega_vec] = ComputeStressesAndSolidBodyMotionOnly(xc_shift,...
                                                                               finalposint_shift,...
                                                                               finalndir,...
                                                                               finalori,...
                                                                               drag_in,...
                                                                               torque_in,...
                                                                               Nf,...
                                                                               background_flow,...
                                                                               M,...
                                                                               flow);

    % Compute the force acting on cube to be able to compute internal stresses
    force_acting_on_cube = ComputeForceActingOnCube(xc,U_vec,Omega_vec,drag_in,torque_in,flow);

    % Get index of internal faces
    internal_faces = FindInternalFaces(xc,NC);
    [center_of_faces, normal_direction, orientation] = BuildBaseCube(finalposint,...
                                                            finalndir,...
                                                            finalori);
    faces_of_base_cube = center_of_faces(1:6,:);
    % Use external stresses to compute internal stresses
    [internal_and_external_stresses,...
        internal_stresses,...
        internal_faces_and_cubes_index_array,...
        internal_faces_and_cubes_index_array_no_double_counting] = ComputeInternalStressesOnly(xc,...
                                                                                               NC,...
                                                                                               stress_outer,...
                                                                                               force_acting_on_cube,...
                                                                                               faces_of_base_cube);       

    [xc_1,xc_2,cubes_that_break] = BreakAggregate(xc,...
        internal_stresses,...
        internal_faces,...
        internal_faces_and_cubes_index_array,...
        internal_faces_and_cubes_index_array_no_double_counting);
    % to be able to compute the distance between where the aggregate break
    % and the center of mass I need to save the internal_faces and cube
    % index array
    % the row whose fourth column = -1 is the breaking face
    % I can pick either the first or the second cube
    % The corresponding face will be located in between the two center of
    % masses (center of each cube)
    % So if I simply pick the halfway-point, that should suffice to give me
    % the location of the face: com(2)-com(1) = face_center
    % finalpos is just for external faces, so it would not work
    
    NC_1 = size(xc_1,1);
    NC_2 = size(xc_2,1);
    % for the new aggregates
    [finalposint_1, finalndir_1, finalori_1,Nf_1] = build_faces(xc_1,size(xc_1,1));
    % Nf is the number of faces
    % finalposint is the location of the center of each face
    % finalndir is the direction of the normal
    % finalori is the orientation of each face
    % Nf is the number of faces

    % for the new aggregates
    [finalposint_2, finalndir_2, finalori_2,Nf_2] = build_faces(xc_2,size(xc_2,1));
    % Nf is the number of faces
    % finalposint is the location of the center of each face
    % finalndir is the direction of the normal
    % finalori is the orientation of each face
    % Nf is the number of faces

    if NC_1==NC
        loop_count = loop_count+1;
        filename_1 = sprintf("~/Desktop/Disaggregation/IAA_23_09_30/Rotation/Stats/SEED_%i_ORIGINAL_Agg_Internal_stresses_NC_%i_FLOW_%i_LoopEvent_%i",SEED,NC,flow,loop_count);     
        save(filename_1,'xc','internal_stresses')
    else
        filename_1 = sprintf("~/Desktop/Disaggregation/IAA_23_09_30/Rotation/Stats/SEED_%i_ORIGINAL_Agg_Internal_stresses_NC_%i_FLOW_%i",SEED,NC,flow);
        save(filename_1,'internal_stresses')
        filename_2 = sprintf("~/Desktop/Disaggregation/IAA_23_09_30/Rotation/Stats/SEED_%i_ORIGINAL_Agg_NC_%i_BROKEN_FLOW_%i",SEED,NC,flow);
        save(filename_2,'xc','xc_1','xc_2','cubes_that_break')
    end
    T_stresses = table(internal_stresses')
    filename = sprintf('~/Desktop/Disaggregation/IAA_23_09_30/Rotation/Stats/Text/SEED_%i_NC_%i_FLOW_%i_Internal_Stresses.txt',SEED,NC,flow)
    writetable(T_stresses,filename,'Delimiter','\t','WriteRowNames',false)
    
    T_xc = table(xc)
    filename = sprintf('~/Desktop/Disaggregation/IAA_23_09_30/Rotation/Stats/Text/SEED_%i_NC_%i_FLOW_%i_xc.txt',SEED,NC,flow)
    writetable(T_xc,filename,'Delimiter','\t','WriteRowNames',false)
    
    T_xc1 = table(xc_1)
    filename = sprintf('~/Desktop/Disaggregation/IAA_23_09_30/Rotation/Stats/Text/SEED_%i_NC_%i_FLOW_%i_xc_1.txt',SEED,NC,flow)
    writetable(T_xc1,filename,'Delimiter','\t','WriteRowNames',false)
    
    T_xc2 = table(xc_2)
    filename = sprintf('~/Desktop/Disaggregation/IAA_23_09_30/Rotation/Stats/Text/SEED_%i_NC_%i_FLOW_%i_xc_2.txt',SEED,NC,flow)
    writetable(T_xc2,filename,'Delimiter','\t','WriteRowNames',false)
    end
end

