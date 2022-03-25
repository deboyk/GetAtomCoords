function [ H, Leverage ] = MolecInfluenceMatrix( M )
% Molecular Influence Matrix
%   This function is used to calculate the molecular influence matrix (M) 
%   from the 3D cartesian coordinates of a molecule. 
%
% See Also: MolecMatrix, R3mCalculate
%
% Author: Kevin DeBoyace
%         Wildfong Lab
%         Duquesne University
% Updated: Jan 2019

H = M*pinv(M'*M)*M';
Leverage = diag(H);
end

