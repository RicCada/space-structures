function Ksim = compTrussK(EAsim)
    syms L
    Ksim = int(EAsim/L^2 * [1 -1; -1 1], 'x', 0, L); 
end

