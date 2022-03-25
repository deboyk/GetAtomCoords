function [ R ] = InfluenceDistanceMat( Atoms, Leverage, G )
%InfluenceDistanceMat (R) - Calculation of the Influence Distance Matrix
%(R)
%   This function calculates the influence/ distance matrix (R). Its initial
%   application was for the calculation of 'R3m', but may be used to
%   calculate any of the R-indices. The function requires input of 'Atoms',
%   'Leverage', and 'G' (the Molecular Geometry Matrix). These values come
%   from the functions 'ImportSDF', 'MolecInfluenceMatrix', and
%   'EuclidDistance', respectively.
%
% See Also: ImportSDF, MolecInfluenceMatrix, EuclidDistance, R3mCalculate
%
% Author: Kevin DeBoyace
%         Wildfong Lab
%         Duquesne University
% Updated: Jan 2019


for ii = 1:size(Atoms,1);
    for jj = ii:size(Atoms,1);
        if jj == ii
            R(ii,jj) =0;
        else
            R(ii,jj) = (sqrt(Leverage(ii)*Leverage(jj))/G(ii,jj));
            R(jj,ii) = (sqrt(Leverage(jj)*Leverage(ii))/G(jj,ii));
        end
    end
end

end

