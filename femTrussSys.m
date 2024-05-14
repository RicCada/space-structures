clear; clc; close all; 

syms L alpha  

LVal = 1000;        % [mm]
EA = 1e6;           % [N]
uVal = 10; 
vVal = -5; 

% compute transformation matrix
T = [cos(alpha) sin(alpha) 0 0; 0 0 cos(alpha) sin(alpha)]; 

% compute stiffness
ki = compTrussK(EA);
K = T.' * ki * T; 
Kfun = matlabFunction(K, 'Vars', [alpha, L]); 


% declare variables and namnes
alphaVec = [0, 90, 180, 45, 135]; 
lVec     = [LVal, LVal, LVal, LVal*sqrt(2), LVal*sqrt(2)]; 
uNames = {'U1', 'V1', 'U2', 'V2', 'U3', 'V3', 'U4', 'V4'}; 
uiNames = { {'U1', 'V1', 'U2', 'V2'}, ... 
            {'U2', 'V2', 'U3', 'V3'}, ... 
            {'U3', 'V3', 'U4', 'V4'}, ...
            {'U1', 'V1', 'U3', 'V3'}, ... 
            {'U2', 'V2', 'U4', 'V4'}};

Fval = [nan, nan, nan, nan; 
        nan, nan, 0, 0; 
        0, 0, nan, nan; 
        nan, nan, 0, 0; 
        nan, nan, nan, nan];

Uval = [0, 0, uVal, vVal, nan, nan, 0, 0]'; 

nT = 5; 
trussArray = cell(1, nT); 

for i = 1:nT
    truss = struct; 

    Kval = Kfun(deg2rad(alphaVec(i)), lVec(i)); 
    Kval(abs(Kval) < 1e-10) = 0; 
   
    truss.K = Kval; 
    truss.F = Fval(i, :)'; 
    truss.uNames = uiNames{i}; 

    trussArray{i} = truss; 
end

[Ktot, Ftot] = syntsStruct(trussArray, uNames); 
[Utot, Ftot] = solveStruct(Ktot, Ftot, Uval);

%% Test 2
clear; clc; close all; 
syms L
LVal = 1000; 
WVal = 1e6; 
PVal = 1000; 
ktVal = 1e8; 
EA = 1e6; 
EJ = 2.5e11; 

Kt = matlabFunction(compTrussK(EA), 'Vars', L); 
Kb = matlabFunction(compBeamK(EJ), 'Vars', L);

truss.K = Kt(LVal); 
truss.uNames = {'V1', 'V2'}; 
truss.F = [nan, 0]'; 

beam1.K = Kb(LVal); 
beam1.uNames = {'V2', 'TH2', 'V3', 'TH3'};
beam1.F = [0, 0, NaN, -WVal]'; 

beam2.K = Kb(LVal); 
beam2.K(4, 4) = beam2.K(4, 4) + ktVal; 
beam2.uNames = {'V3', 'TH3', 'V4', 'TH4'}; 
beam2.F = [NaN, 0, PVal, 0]'; 

els = {truss, beam1, beam2}; 

uNames = {'V1', 'V2', 'TH2', 'V3', 'TH3', 'V4', 'TH4'}; 

Uval = [0, NaN, NaN, 0, NaN, NaN, NaN]'; 

[Ktot, Ftot] = syntsStruct(els, uNames); 
[Ures, Fres] = solveStruct(Ktot, Ftot, Uval);





