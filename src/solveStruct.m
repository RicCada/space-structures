function [U, F] = solveStruct(K, F0, U0)

    idxDispl = isnan(F0);
    idxForce = not(idxDispl);  

    Kred = K(idxForce, idxForce); 
    Fred = F0(idxForce) - K(idxForce, idxDispl) * U0(idxDispl) ;
    
    Ured = Kred\Fred; 
    
    U = U0; 
    U(idxForce) = Ured; 

    F = K * U; 
end

