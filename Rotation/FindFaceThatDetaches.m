function [face_detached] = FindFaceThatDetaches(internal_stresses,...
                                                indices_of_internal_stresses)
                                                                  
% for now we assume that the aggregate will break only in two parts
% based on which face is subject to the max internal stress
max_stress = 0;
for nn = 1:size(internal_stresses,2)
    max_found = max(norm(internal_stresses(:,nn)));
    if max_found>max_stress
        max_stress = max_found;
        index = nn;
        face_detached = indices_of_internal_stresses(index);
    end
end

end
