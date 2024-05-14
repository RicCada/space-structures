function [Ksim] = compBeamK(EJsim)
    % K sim consider displacements in the following order: 
    % [w1, theta1, w2, theta2]


    syms w1 w2 L th1 th2 x; 

    X = [1 0 0 0; 1 L L^2 L^3; 0 1 0 0; 0 1 2*L 3*L^2]; 
    U = [w1 w2 -th1 -th2].'; 
    
    cVect = X\U; 
    
    % reconstruct w0 expression
    xVect = [1, x, x^2, x^3]; 
    wFun = xVect * cVect; 
    
    % consider: 
    % w = N1*w1 + N2*th1 + N3*w2 + N4*th2
    % compute N1, N2, N3, N4
    
    [Nvec, T] = coeffs(wFun, [w1, th1, w2, th2]); 
    
    % derive the coresponding stiffenss matrix
    % K = int(0, L) Nvec/xx * EJ * Nvec/xx
    
    NvecDD = diff(diff(Nvec, 'x'), 'x'); 
    
    Ksim = int(NvecDD.' * NvecDD * EJsim, 'x', 0, L);
end

