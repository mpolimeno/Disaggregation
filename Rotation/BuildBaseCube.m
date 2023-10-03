function [center_of_faces, normal_direction, orientation] = BuildBaseCube(center_of_faces,...
                                                            normal_direction,...
                                                            orientation)
%{
                                                            
Function builds base cube to give faces proper orientations and find internal faces moving forward.
Improves modularity
                                                            
center_of_face holds the position of the center of each face
                                                            
normal_direction holds the corresponding normal direction (x,y or z -> 1,2,3)
                                                            
orientation gives whether it is inward (-1) or outward (1)
                                                            
%}
                                                        
center_of_faces(1,:) = [1,0,0];
center_of_faces(2,:) = [0,1,0];
center_of_faces(3,:) = [0,0,1];
center_of_faces(4,:) = [0,0,-1];
center_of_faces(5,:) = [0,-1,0];
center_of_faces(6,:) = [-1,0,0];

normal_direction(1) = 1;
orientation(1)      = 1;

normal_direction(2) = 2;
orientation(2)      = 1;

normal_direction(3) = 3;
orientation(3)      = 1;

normal_direction(4) = 3;
orientation(4)      = -1;

normal_direction(5) = 2;
orientation(5)      = -1;

normal_direction(6) = 1;
orientation(6)      = -1;

end