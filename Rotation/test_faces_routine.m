% testing some routines
close all
clear
clc

% impose flow
% flow = 1 -> Translation
% flow = 2 -> Rotation
% flow = 3 -> Extension
% flow = 4 -> Shear
flow = 4;

method = 1;
if method==1 % IAA routine
    SEED = 1;
    number_of_cubes_in_aggregate = 200;
    center_of_cubes_in_aggregate = DLA_3D(number_of_cubes_in_aggregate,SEED);
elseif method==2 % Build dumbell whose bar lies on y-axis: useful to test rotational flow around vertical
    number_of_cubes_in_aggregate_dumb = 4;
    number_of_cubes_in_aggregate_bell = 54;
    pos_dumb = build_dumb(number_of_cubes_in_aggregate_dumb,1);
    [pos_bell_l,pos_bell_r] = build_bell(number_of_cubes_in_aggregate_bell,1);
    center_of_cubes_in_aggregate = [pos_bell_l;pos_dumb;pos_bell_r];
    number_of_cubes_in_aggregate = size(center_of_cubes_in_aggregate,1);
elseif method==3 % build a vertical bar: useful for shear case u=gamma*y*i_hat
    center_of_cubes_in_aggregate = [0,0,0;...
                                    0,0,2;...
                                    0,0,-2;...
                                    0,0,-4;...
                                    0,0,4;...
                                    0,0,-6];
    number_of_cubes_in_aggregate = size(center_of_cubes_in_aggregate,1);
else
    msg = "'method' needs to be either 1, 2 or 3";
    error(msg);
end

if size(center_of_cubes_in_aggregate,1)==1
    center_of_mass_of_aggregate = center_of_cubes_in_aggregate; % for one cube only
else
    center_of_mass_of_aggregate = mean(center_of_cubes_in_aggregate); % for multiple cubes
end

% unshifted routine
[center_of_external_faces, normal_direction_of_each_external_face, orientation_of_each_external_face,number_of_external_faces] = build_faces(center_of_cubes_in_aggregate,number_of_cubes_in_aggregate);

[center_of_faces, normal_direction, orientation] = BuildBaseCube(center_of_external_faces,...
                                                            normal_direction_of_each_external_face,...
                                                            orientation_of_each_external_face);
% Plot the external faces of the original aggregate
plot_faces(center_of_external_faces,...
           normal_direction_of_each_external_face,...
           number_of_external_faces,1,[0 0.4470 0.7410]);

if flow==1
    extensional_matrix     = [];
    velocity_vector_input  = [0;0;0];
    translational_velocity = [0;0;1];
    [forceout,drag,torque] = fractal_bi_stokes_force(center_of_cubes_in_aggregate,...
                                                     center_of_external_faces,...
                                                     normal_direction_of_each_external_face,...
                                                     orientation_of_each_external_face,...
                                                     translational_velocity,...
                                                     number_of_external_faces);
    drag_imposed   = drag;
    torque_imposed = [0;0;0];
elseif flow==2
    extensional_matrix     = [];
    rotational_velocity    = [0;0;1];
    velocity_vector_input  = [0;0;0];
    [forceout,drag,torque] = fractal_bi_stokes_force_rot(center_of_cubes_in_aggregate,...
                                                         center_of_external_faces,...
                                                         normal_direction_of_each_external_face,...
                                                         orientation_of_each_external_face,...
                                                         rotational_velocity,...
                                                         number_of_external_faces);
    drag_imposed = [0;0;0];
    torque_imposed = torque;
elseif flow==3
   extensional_matrix    = [-1,0,0;0,1,0;0,0,0];
   velocity_vector_input = [0;0;0];
   drag_imposed          = [0;0;0];
   torque_imposed        = [0;0;0];
elseif flow==4
   % this shear flow is u = gamma*y*i_hat
   shear_rate            = 1; % gamma_t in the notes
   extensional_matrix    = 0.5*[0,shear_rate,0;shear_rate,0,0;0,0,0];
   k_hat                 = [0;0;1];
   velocity_vector_input = -0.5*shear_rate*k_hat; % this is the omega for strain
   drag_imposed          = [0;0;0];
   torque_imposed        = [0;0;0];
end

internal_faces = FindInternalFaces(center_of_cubes_in_aggregate,number_of_cubes_in_aggregate);
external_faces = (6*number_of_cubes_in_aggregate) - size(internal_faces,2);

center_of_external_faces_shift = center_of_external_faces - center_of_mass_of_aggregate;
center_of_cubes_in_aggregate_shift = center_of_cubes_in_aggregate - center_of_mass_of_aggregate; % shift for stresses
[LHS,sol,stress_outer,U_vec,Omega_vec] = ComputeStressesAndSolidBodyMotionOnly(center_of_cubes_in_aggregate_shift,...
                                                                                center_of_external_faces_shift,...
                                                                                normal_direction_of_each_external_face,...
                                                                                orientation_of_each_external_face,...
                                                                                drag_imposed,...
                                                                                torque_imposed,...
                                                                                number_of_external_faces,...
                                                                                velocity_vector_input,...
                                                                                extensional_matrix,...
                                                                                flow);
% shift back for appropriate count of internal faces
% the space location is not needed to compute these stresses
% we just neeed to know the index of the internal faces

% get the force acting on each cube
force_cube = ComputeForceActingOnCube(center_of_cubes_in_aggregate,drag_imposed,flow);
faces_of_base_cube = center_of_faces(1:6,:);
[internal_and_external_stresses,...
    internal_stresses,...
    internal_faces_complete,...
    internal_faces_unrepeated] = ComputeInternalStressesOnly(center_of_cubes_in_aggregate,...
                                                             number_of_cubes_in_aggregate,...
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
[center_of_cubes_in_aggregate_1,center_of_cubes_in_aggregate_2,cubes_that_break] = BreakAggregate(center_of_cubes_in_aggregate,...
                                                                                                  internal_stresses,...
                                                                                                  internal_faces,...
                                                                                                  internal_faces_complete,...
                                                                                                  internal_faces_unrepeated);

[internal_faces_pos,normal_direction_of_each_external_face_internal,orientation_of_each_external_face_internal] = GetPositionOfInternalFaces(center_of_cubes_in_aggregate, number_of_cubes_in_aggregate);

position_directions_faces_and_stresses = zeros(size(normalized_stresses,1),6);
for ii=1:size(normalized_stresses,1)
    position_directions_faces_and_stresses(ii,1:3) = internal_faces_pos(ii,:);
    position_directions_faces_and_stresses(ii,4) = normal_direction_of_each_external_face_internal(ii);
    position_directions_faces_and_stresses(ii,5) = internal_faces(ii);
    position_directions_faces_and_stresses(ii,6) = normalized_stresses(ii);
end
input_array = sortrows(position_directions_faces_and_stresses,6);

% Plot the internal stresses of the original aggregate
PlotInternalStresses(2,...
                     center_of_external_faces,...
                     normal_direction_of_each_external_face,...
                     stress_and_faces,...
                     input_array)

% Build the external faces of aggregate 1
[center_of_external_faces_1, normal_direction_of_each_external_face_1, orientation_of_each_external_face_1,number_of_external_faces_1] = build_faces(center_of_cubes_in_aggregate_1, size(center_of_cubes_in_aggregate_1,1));
% Plot the external faces of aggregate 1
plot_faces(center_of_external_faces_1, normal_direction_of_each_external_face_1,number_of_external_faces_1,3,[0.3010 0.7450 0.9330]);

% Build the external faces of aggregate 2
[center_of_external_faces_2, normal_direction_of_each_external_face_2, orientation_of_each_external_face_2,number_of_external_faces_2] = build_faces(center_of_cubes_in_aggregate_2, size(center_of_cubes_in_aggregate_2,1));
% Plot the external faces of aggregate 2
plot_faces(center_of_external_faces_2, normal_direction_of_each_external_face_2,number_of_external_faces_2,4,[0.4660 0.6740 0.1880]);

% Will have to move it to top and make it more modular
% Trying to get a hang of the streamlines
L_x = max(abs(center_of_external_faces(:,1)));
L_y = max(abs(center_of_external_faces(:,2)));
L_z = max(abs(center_of_external_faces(:,3)));
[X,Y,Z] = meshgrid(-L_x:L_x,-L_y:L_y,-L_z:L_z);
%[start_x,start_y,start_z] = meshgrid(-5:5,-5:5,-5:5);
U = shear_rate*Y;   % this should be it as u=gamma*y*i_hat
V = zeros(size(Y)); % 0
W = zeros(size(Y)); % 0

% plot the maximum stress only
PlotMaximumInternalStress(5,...
                          center_of_external_faces,...
                          normal_direction_of_each_external_face,...
                          stress_and_faces,...
                          input_array,...
                          face_detached)
hold on
quiver3(X,Y,Z,U,V,W)
% verts = stream3(X,Y,Z,U,V,W,start_x,start_y,start_z);
% lineobj = streamline(verts);
view(3)
xlim([-L_x,L_x])
ylim([-L_y,L_y])
zlim([-L_z,L_z])