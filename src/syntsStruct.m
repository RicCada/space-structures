function [Ktot, Ftot] = syntsStruct(elements, uNames)
    
    nE = length(elements);
    nT = length(uNames);        % size of the finla system
    
    % init matrices
    Ktot = zeros(nT); 
    Ftot = zeros(nT, 1); 
    
    % Build stiffness matrix
    for i = 1:nT
        rN = uNames{i};         % row name
        for j = 1:nT
            cN = uNames{j};         % col name

            % scan each single element
            for k = 1:nE
                elStruct = elements{k}; 
                idRow = findPosition(elStruct.uNames, rN); 
                idCol = findPosition(elStruct.uNames, cN); 
                
                if idRow*idCol ~= 0
                    Ktot(i, j) = Ktot(i, j) + elStruct.K(idRow, idCol); 
                end
            end
        end
    end

    % Build force matrix
    for i = 1:nT
        rN = uNames{i}; 
        for k = 1:nE
            elStruct = elements{k}; 
            idRow = findPosition(elStruct.uNames, rN); 
            if idRow ~= 0
                Ftot(i) = Ftot(i) + elStruct.F(idRow); 
            end
        end
    end
    

end


function id = findPosition(list, target)
    n = length(list);
    id = 0; 

    for i = 1:n
        if isequal(list{i}, target)
            id = i; 
            return; 
        end
    end
end