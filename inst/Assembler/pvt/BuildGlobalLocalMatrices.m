function [id, lm] = ...
    BuildGlobalLocalMatrices(element_local_mapping,NUMBER_OF_SPACE_D)
    [nen, nel] = size(element_local_mapping);
    MAX_BASIS = max(max(element_local_mapping));
    id = zeros(MAX_BASIS,NUMBER_OF_SPACE_D);
    for i=1:NUMBER_OF_SPACE_D
        id(:,i) = [(i-1)*(MAX_BASIS)+1:(i)*MAX_BASIS]';
    end
    lm = zeros(NUMBER_OF_SPACE_D*nen,nel);
    for i=1:nel
        lm(:,i) = reshape(id(element_local_mapping(:,i),:), ...
            NUMBER_OF_SPACE_D*nen,1);
    end
end