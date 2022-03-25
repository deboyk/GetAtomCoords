function [ weightedmass, MolecWeight, colorset ] = AtomicWeighting( Atoms )
%AtomicWeigthing A function to normalize atomic masses to carbon
%   This function is used to determine the appropriate atomic mass
%   weighting for the calculation of R3m. In the calculation, the atomic
%   mass is normalized to carbon.
%   This function also assigns a color to each atom for plotting (see
%   PlotMolecule) and calculates the molecular weight. 
%   
%   Author: Kevin DeBoyace
%   Updated: Jan 2019
%
%   See also: R3mCalculate, PlotMolecule

atomicweight = zeros(size(Atoms,1),1);

% Check if molecule contains Hydrogens. If not, display warning message.
ContainHydrogen = strcmp(Atoms, 'H');
if any(ContainHydrogen) ~= 1;
    h = warndlg('WARNING: This molecule is missing Hydrogens. Hydrogens should be added to accurately calculate molecular descriptors.', 'WARNING');
    waitfor(h) %Pause code until user closes warning message.
else
end
% Calculate molecular weight and select atom colors for those Atoms which
% are present.
for b = 1:(size(Atoms,1));
    if strcmp(Atoms(b,:), 'O ') || strcmp(Atoms(b,:), 'O') == 1;
        atomicweight(b) = 15.9994;
        colorset(b,:) = [1 0 0];
    elseif strcmp(Atoms(b,:),'C ')|| strcmp(Atoms(b,:), 'C') == 1;
        atomicweight(b) = 12.0107;
        colorset(b,:) = [0 0 0];
    elseif strcmp(Atoms(b,:),'N ')|| strcmp(Atoms(b,:), 'N') == 1;
        atomicweight(b) = 14.0067;
        colorset(b,:) = [0 0 1];
    elseif strcmp(Atoms(b,:),'Cl') || strcmp(Atoms(b,:),'Cl ') == 1;
        atomicweight(b) = 35.453;
        colorset(b,:) = [0 0.3 0];
    elseif strcmp(Atoms(b,:),'S ')|| strcmp(Atoms(b,:), 'S') == 1;
        %atomicweight(b) = 32.065;
        atomicweight(b) = 32.066;
        colorset(b,:) = [1 1 0];
    elseif strcmp(Atoms(b,:), 'P ') || strcmp(Atoms(b,:), 'P') == 1;
        atomicweight(b) = 30.973762;
        colorset(b,:) = [1 0.5 0];
    elseif strcmp(Atoms(b,:), 'F ') || strcmp(Atoms(b,:), 'F') == 1;
        atomicweight(b) = 18.998;
        colorset(b,:) = [0.5 0.8 0.4];
    elseif strcmp(Atoms(b,:), 'H ') || strcmp(Atoms(b,:), 'H') == 1;
        atomicweight(b) = 1.00794;
        colorset(b,:) = [0.9 0.9 0.9];
    end
end

% Throw an error if an atom is not listed in the code above. If this error
% is thrown, the atom needs to be added to the portion of the code above to
% ensure correct calculation of Molecular weight and weightedmass. 
for b = 1:(size(Atoms,1));
    if atomicweight(b) == 0;
        error('Error: Atom missing from list in code')
    end
end

atomicweight = atomicweight';
MolecWeight = sum(atomicweight);

% Masses are weighted with respect to the carbon atom : mass/ mass carbon
weightedmass = atomicweight / 12.0107;
weightedmass = weightedmass';

end

