function [R3m] = R3mCalculate_auto(filename)
%R3mCalculate - A function to calculate the R3m molecular descriptor. 
%   This function calculates the molecular matrix (M), geometry matrix (G),
%   molecular influence matrix (H), and influence distance matrix (R). This
%   is accomplished by calling the following subfunctions: ImportSDF,
%   MolecMatrix, EuclidDistance, MolecInfluenceMatrix, and
%   InfluenceDistanceMat. A GUI is called to select the file type to import
%   (filetypegui). Next, the molecule is plotted using the PlotMolecule
%   function to ensure the 3D coordinates are reasonable. Finally, the R3m
%   value is calculated.
%
%   See also: IMPORTSDF, MOLECMATRIX, EUCLIDDISTANCE, MOLECINFLUENCEMATRIX,
%   INFLUENCEDISTANCEMAT, PLOTMOLECULE, ATOMICWEIGHTING, ATOMCONNECTION
%
%   References: 
%       [1] Todeschini R, Consonni V. 2008. Handbook of molecular 
%           descriptors. ed.: John Wiley & Sons.
%       [2] Consonni V, Todeschini R, Pavan M, Gramatica P 2002. 
%           Structure/response correlations and similarity/diversity 
%           analysis by GETAWAY descriptors. 2. Application of the novel 
%           3D molecular descriptors to QSAR/QSPR studies. Journal of 
%           chemical information and computer sciences  42(3):693-705.
%
%
%   Code Author: Kevin DeBoyace
%                Duquesne University
%   
%   Date:   October 2018

%% Matrices
% M - Molecular Matrix
% G - Geometry Matrix
% H - Molecular Influence Matrix
% R - Influence Distance Matrix

% Molecular Matrix
[ M, Atoms, numatoms, Connectivity] = MolecMatrix( filename );

% Atomic weighting
[ weightedmass, MolecWeight, colorset ] = AtomicWeighting( Atoms );
% Geometry Matrix
[G] = EuclidDistance(M);
% Molecular Influence Matrix
[ H, Leverage ] = MolecInfluenceMatrix( M );
% Influence Distance Matrix
[ R ] = InfluenceDistanceMat( Atoms, Leverage, G );

%% Atom connections
% Find atoms 1, 2 and 3 bond distances away.
[ connections, topology_mat ] = AtomConnection( Connectivity, Atoms );

%% R3m calculation
%  R3m = sum(i to A-1)sum(j>i) R*wi*wj*delta(k;dij) k =1,2..., d
R3m = 0;
for ii = 1:size(Atoms,1)-1;
    for jj = ii+1:size(Atoms,1);
        if topology_mat(ii,jj) == 1;
            calc = R(ii,jj)*weightedmass(ii)*weightedmass(jj);
            R3m = calc + R3m;
        else
        end
    end
end

% R3m

end
