function [xc_1,xc_2,cubes_that_break] = BreakAggregate(xc,...
    internal_stresses,...
    indices_of_internal_stresses,...
    internal_faces_complete,...
    internal_faces_unrepeated)
    %{
        xc: {Number_of_cubes x 3 array} of integers. Each row contains the center of each cube
        in the aggregate
        internal_stresses: {array of floats}. Its size is
        Number_of_internal_faces x 3 (stresses are vectors)
        indices_of_internal_stresses: {array}. Gives the index of the internal
        faces, effectively
        internal_faces_complete: {Number_of_internal_faces by 3
        array}:
        [face number, first cube, second cube]
        internal_faces_unrepeated: {array}. Same
        idea as its almost-namesake, but we avoid double counting the attached
        faces
    %}

    % This function returns the face where the maximum internal stress is found
    [face_detached] = FindFaceThatDetaches(internal_stresses,indices_of_internal_stresses);

    [cubes_that_break,where] = GetCubesWhoseFaceDetaches(face_detached,...
        internal_faces_complete,...
        internal_faces_unrepeated);


    % Matlab has two very helpful functions
    % Build graph
    s = internal_faces_unrepeated(:,2);
    t = internal_faces_unrepeated(:,3);
    G = graph(s,t);

    % Break Graph at the edge found to have the maximum stress
    S = rmedge(G,where);
    % Find Connected components
    bins = conncomp(S);
    % Store nodes into cell array
    binnodes = accumarray(bins',1:numel(bins),[],@(v) {sort(v')});

    % The first array will be Aggregate 1
    xc_1 = xc(binnodes{1}',:);
    % Handle looped cases
    if size(binnodes,1)>1
        % The second array will be Aggregate 2
        xc_2 = xc(binnodes{2}',:);
    else
        xc_2 = []; % If looped, do not break
    end

end
