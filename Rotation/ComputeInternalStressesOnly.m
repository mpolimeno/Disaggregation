function [internal_and_external_stresses,...
    internal_stresses,...
    internal_faces_and_cubes_index_array,...
    internal_faces_and_cubes_index_array_no_double_counting] = ComputeInternalStressesOnly(xc,...
                                                                                            NC,...
                                                                                            forceout,...
                                                                                            force_cube,...
                                                                                            posintb)
    % This function computes the stresses acting on every face of a cluster, even the inner forces
    % it takes as inputs the position of the blocks of the cluster, and the outer forces.

    % Our system of equations making up the constraint will have 6*number_of_cubes unknowns, every face
    % For every outside face, we have the equation is F = forceout given
    % For every internal force, we have force on one side is minus force on the other
    % For every block, we have the equation: sum of forces on block is (sum of all forces)/NC
    % And then, we want to minimize the sum of the square of the forces (if there are degrees of freedom left)

    % The columns will be given as x y z for every face of each block,
    % then faces centered at [1,0,0], [0,1,0], [0,0,1], [0,0,-1], [0,-1,0], [-1,0,0], for every block

    % That will lead to a lot of equations, some unnecessarily, maybe we can avoid that easily?
    % In general, to solve: Min ||Ax - b||^2 subject to Cx = d, one solves:
    % https://stanford.edu/class/ee103/lectures/constrained-least-squares_slides.pdf
    % [Matteo] http://www.seas.ucla.edu/~vandenbe/133A/lectures/cls.pdf
    % [Matteo] https://pdfs.semanticscholar.org/2ad9/a8e86735e8104f09fe15f02584ddf8a344dd.pdf
    % [ 2 A^T A    C^T ] [x] = [2 A^T b ]
    % [   C        0   ] [z] = [ d    ]
    % where z are Lagrange multipliers
    % In our case, to minimize the square of the forces, A = I and b = 0 so we have
    % [ 2 I    C^T ] [x] = [ 0 ]
    % [  C       0 ] [z] = [ d ]
    % So we need to define the Matrix C and the RHS d .
    
    % outers is an array of flags:
    % If it is 2, it is an outer face,
    % if it is 1, it is an inner face already dealt with
    % if it is 0, it is a new inner face
    % flags for each face
    outers = zeros(1,NC*6); % to check if we have an outer face or an inner face
    
    indu = 0;
    inde = 0;
    jnf = 0; % index for outer faces
    
    internal_faces_and_cubes_index_array = [];
    internal_faces_and_cubes_index_array_no_double_counting = [];
    for i=1:NC
        % next 18 unknowns have to do with block i
        if (i<NC)
            % Add equation matching total forces on each block to average force per block
            inde = inde + 3;
                %jnf = jnf + 1; % this is the index of the outer faces
            for jj=1:3 % x, y, z
                for kk=1:6 % go over everyface
                    Cmat(inde-3+jj,indu+3*(kk-1)+jj) = 1;
                end
                dRHS(inde-2:inde) = force_cube(i,:)/4; % 4 is the area of a square face
            end
        end
        % now add matching equation for each face
        for k = 1:6
            indu = indu + 3; % update index of faces, indexes of unknowns
            % even for the periodic case now the neighbors need to be adjusted,
            % as they could be a period away
            neibpos = xc(i,:)+2*posintb(k,:);
            if (outers(indu/3) == 0) % never been to that face
                for j=1:NC % check if it is an outer face
                    if ( neibpos == xc(j,:)) % We have a neighbor
                        outers(indu/3) = 1; % this face was visited and is an inner face
                        internal_faces_and_cubes_index_array_no_double_counting = [internal_faces_and_cubes_index_array_no_double_counting; (indu/3) i j];
                        internal_faces_and_cubes_index_array = [internal_faces_and_cubes_index_array; (indu/3) i j];
                        indu2 = (6*(j-1)+7-k)*3; % index of the corresponding inner face we just found
                        outers(indu2/3) = 1; % the corresponding face was visited too
                        internal_faces_and_cubes_index_array = [internal_faces_and_cubes_index_array; (indu2/3) i j];
                        inde = inde + 3; % we get 3 equations from an inner face
                        Cmat(inde-2,indu-2) = 1;
                        Cmat(inde-1,indu-1) = 1;
                        Cmat(inde,indu) = 1;
                        Cmat(inde-2,indu2-2) = 1;
                        Cmat(inde-1,indu2-1) = 1;
                        Cmat(inde,indu2) = 1;
                        dRHS(inde-2:inde) = 0;
                        break; % because we found our neighbor, we can stop looking not sure if that works ??
                    end
                end
            if (outers(indu/3)==0) % We have no neighbor
                outers(indu/3) = 2; % this face is an outer face
                inde = inde + 3; % we get 3 equations from an outer face
                Cmat(inde-2,indu-2) = 1;
                Cmat(inde-1,indu-1) = 1;
                Cmat(inde,indu) = 1;
                jnf = jnf + 1;
                dRHS(inde-2:inde) = forceout(:,jnf);
            end
            end

        end
    end

    % outers is an array of flags -> if it stores 1, that should be an internal
    % face -> test the case for 2 cubes


    % build the matrix to minimize internal forces
    sC = size(Cmat);
    % make a note on minimization
    Amat = [2*eye(sC(2)) transpose(Cmat); Cmat zeros(sC(1),sC(1))];
    bmat = [zeros(sC(2),1); dRHS'];
    size(dRHS)
    tic
    r_amat = rank(Amat,eps);
    if r_amat==size(Amat,2)
        xmat = Amat\bmat; %Matteo: These are the external stresses and velocities;
    else
        xmat = pinv(Amat)*bmat;%lsqminnorm(Amat,bmat,eps);%
    end
    toc

    endtime_internal = toc;
    disp(endtime_internal)

    size(xmat)

    internal_and_external_stressestemp = xmat(1:sC(2));

    %faces = NC*6;
    %internal_and_external_stresses = zeros(faces,3);
    for i=1:NC*6
        internal_and_external_stresses(:,i) = internal_and_external_stressestemp(3*(i-1)+1:3*i);
    end

    internal_faces = FindInternalFaces(xc,NC);

    internal_stresses = internal_and_external_stresses(:,internal_faces);
end
