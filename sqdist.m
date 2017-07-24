function [ Dist ] = sqdist( A , B )
% sqdist: Squared Euclidean distances between columns of A and B.
% Thus, data points are assumed to be in columns, not rows
%
% Dist = sqdist( A, B )
% Dist = sqdist( A ) assumes B = A

% Farhad Pourkamali-Anaraki, E-mail: Farhad.Pourkamali@colorado.edu
% University of Colorado Boulder

%{
Inputs: 
    - A: input matrix of size p*n1
    - B: input matrix of size p*n2 
Outputs: 
    - Dist: matrix of pairwise distances of size n1*n2
%}

if nargin < 2
    % Assume B = A
    AA = sum(A.^2,1); 

    Dist = -2*(A'*A);
    Dist = bsxfun( @plus, Dist, AA' );
    Dist = bsxfun( @plus, Dist, AA );

else
    
    [p1,~] = size(A);
    [p2,~] = size(B);
    if p1~=p2, error('A and B should have the same number of rows'); end
    
    AA = sum(A.^2,1);
    BB = sum(B.^2,1);

    
    Dist = -2*(A'*B);
    Dist = bsxfun( @plus, Dist, AA' );
    Dist = bsxfun( @plus, Dist, BB );
    % fprintf('Discrepancy between the two approaces: %.2e\n', norm(Dist-Dist_Ref,'fro')/norm(Dist_Ref,'fro') );
end
end

