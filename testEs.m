clear; clc; close all; 

syms alpha L


EA = 210e3 * 300;  % [MPa * mm2]
P = 10e3;          % [N]

lVals = [412.31, 565.68, 300, 412.31]; 
alphaVals = [14, 225, 270, 346]*pi/180; 

uMat = { {'U1', 'V1', 'U2', 'V2'}; ...
         {'U2', 'V2', 'U3', 'V3'}; ...
         {'U1', 'V1', 'U3', 'V3'}; ...
         {'U3', 'V3', 'U4', 'V4'}}; 

FMat = [P, 0, nan, nan; 
        nan, nan, 0, 0; 
        0, 0, 0, 0; 
        0, 0, nan, nan]; 

U0 = [nan, nan, 0, 0, nan, nan, 0, 0]'; 

uNames = {'U1', 'V1', 'U2', 'V2', 'U3', 'V3', 'U4', 'V4'}; 


T = [cos(alpha), sin(alpha), 0, 0; 0, 0, cos(alpha), sin(alpha)]; 

k = compTrussK(EA); 
Ki = T.' * k * T; 
KiFun = matlabFunction(Ki, 'Vars', [L, alpha]); 

elements = cell(1, 4); 

for i = 1:4
    el = struct; 
    el.K = KiFun(lVals(i), alphaVals(i)); 
    el.F = FMat(i, :)'; 
    el.uNames = uMat{i}; 

    elements{i} = el; 
end



[Ktot, Ftot] = syntsStruct(elements, uNames); 
[Ures, Fres] = solveStruct(Ktot, Ftot, U0);