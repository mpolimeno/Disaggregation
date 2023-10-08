% testing some routines
close all
clear
clc

method = 1;
if method==1
    SEED = 1;
    NC = 100;
    xc = DLA_3D(NC,SEED);
elseif method==2
    NC_dumb = 4;
    NC_bell = 54;
    pos_dumb = build_dumb(NC_dumb,1);
    [pos_bell_l,pos_bell_r] = build_bell(NC_bell,1);
    xc = [pos_bell_l;pos_dumb;pos_bell_r];
    NC = size(xc,1);
else
    msg = "'method' needs to be either 1 or 2";
    error(msg);
end

if size(xc,1)==1
    cm = xc; % for one cube only
else
    cm = mean(xc); % for multiple cubes
end

% unshifted routine
[finalposint, finalndir, finalori,Nf] = build_faces(xc,NC);

[center_of_faces, normal_direction, orientation] = BuildBaseCube(finalposint,...
                                                            finalndir,...
                                                            finalori);

plot_faces(finalposint, finalndir, finalori,Nf,1,'c');

% impose flow
flow = 2;
if flow==3
   M = [-1,0,0;0,1,0;0,0,0];
   u_infty = [0;0;0];
   drag_in = [0;0;0];
   torque_in = [0;0;0];
elseif flow==4
   M = [0,1,0;1,0,0;0,0,0];
   u_infty = 0.5*[0;0;-1];
   drag_in = [0;0;0];
   torque_in = [0;0;0];
elseif flow==1
    M = [];
    u_infty = [0;0;0];
    Uvec = [0;0;1];
    [forceout,drag,torque] = fractal_bi_stokes_force(xc,finalposint,finalndir,finalori,Uvec,Nf);
    drag_in = drag;
    torque_in = [0;0;0];    
elseif flow==2
    M = [];
    Uvec = [0;0;1];
    Rot = Uvec;
    u_infty = [0;0;0];
    [forceout,drag,torque] = fractal_bi_stokes_force_rot(xc,finalposint,finalndir,finalori,Rot,Nf);
    drag_in = [0;0;0];
    torque_in = torque;
end

internal_faces = FindInternalFaces(xc,NC);
external_faces = (6*NC) - size(internal_faces,2);

finalposint_shift = finalposint - cm;
xc_shift = xc - cm; % shift for stresses
[LHS,sol,stress_outer,U_vec,Omega_vec] = ComputeStressesAndSolidBodyMotionOnly(xc_shift,...
                                                                                finalposint_shift,...
                                                                                finalndir,...
                                                                                finalori,...
                                                                                drag_in,...
                                                                                torque_in,...
                                                                                Nf,...
                                                                                u_infty,...
                                                                                M,...
                                                                                flow);
% shift back for appropriate count of internal faces
% the space location is not needed to compute these stresses
% we just neeed to know the index of the internal faces

% get the force acting on each cube
force_cube = ComputeForceActingOnCube(xc,drag_in,flow);
faces_of_base_cube = center_of_faces(1:6,:);
[internal_and_external_stresses,...
    internal_stresses,...
    internal_faces_complete,...
    internal_faces_unrepeated] = ComputeInternalStressesOnly(xc,...
                                                             NC,...
                                                             stress_outer,...
                                                             force_cube,...
                                                             faces_of_base_cube);
    
norm_of_internal_stresses = zeros(size(internal_stresses,2),1);
for stress=1:size(internal_stresses,2)
    norm_of_internal_stresses(stress) = norm(internal_stresses(:,stress));
end

maximum_internal_stress = max(norm_of_internal_stresses);
normalized_stresses = norm_of_internal_stresses/maximum_internal_stress;


stress_and_faces = zeros(size(normalized_stresses,1),2);
for ii=1:size(normalized_stresses,1)
    stress_and_faces(ii,1) = normalized_stresses(ii);
    stress_and_faces(ii,2) = internal_faces(ii);
end
stress_and_faces = sortrows(stress_and_faces);


[face_detached] = FindFaceThatDetaches(internal_stresses,internal_faces);
[xc_1,xc_2,cubes_that_break] = BreakAggregate(xc,...
    internal_stresses,...
    internal_faces,...
    internal_faces_complete,...
    internal_faces_unrepeated);

[internal_faces_pos,finalndir_internal,finalori_internal] = GetPositionOfInternalFaces(xc, NC);

position_directions_faces_and_stresses = zeros(size(normalized_stresses,1),6);
for ii=1:size(normalized_stresses,1)
    position_directions_faces_and_stresses(ii,1:3) = internal_faces_pos(ii,:);
    position_directions_faces_and_stresses(ii,4) = finalndir_internal(ii);
    position_directions_faces_and_stresses(ii,5) = internal_faces(ii);
    position_directions_faces_and_stresses(ii,6) = normalized_stresses(ii);
end
input_array = sortrows(position_directions_faces_and_stresses,6);
PlotInternalStresses(3,...
                     finalposint,...
                     finalndir,...
                     stress_and_faces,...
                     input_array)

[finalposint_1, finalndir_1, finalori_1,Nf_1] = build_faces(xc_1, size(xc_1,1));
plot_faces(finalposint_1, finalndir_1, finalori_1,Nf_1,6,'b');

[finalposint_2, finalndir_2, finalori_2,Nf_2] = build_faces(xc_2, size(xc_2,1));

plot_faces(finalposint_2, finalndir_2, finalori_2,Nf_2,7,'r');

