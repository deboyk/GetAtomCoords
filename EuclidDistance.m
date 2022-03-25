function [G] = EuclidDistance(M)
%G - Determine Euclidean Distance of atoms in a molecule
%   G = Geometric Distance Matrix
%   This function determines the Euclidean distance between all of the
%   atoms in the molecule and outputs the data to the matrix G, which is
%   known as the geometric distance matrix. 
%   See also: R3mCalculate, CheckCentroid, 
%
% Author: Kevin DeBoyace
%         Wildfong Lab
%         Duquesne University
% Updated: Jan 2019


%G = pdist2(M,M,'euclidean'); % Performs same as below, but is slower

for ii = 1:size(M,1);
    for jj = 1:size(M,1);
        % Calculate Euclidean Distance
        G(ii,jj) = sqrt((M(ii,1)-M(jj,1)).^2 + (M(ii,2)-M(jj,2)).^2 +(M(ii,3)-M(jj,3)).^2);
    end
end

end