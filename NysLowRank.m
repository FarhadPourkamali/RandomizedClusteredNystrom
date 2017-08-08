function [ L ] = NysLowRank( X , Z , desiredRank , kernel )
% NysLowRank: Low-rank approximation using the landmark set Z. 
%
% Warning: Please use NysDecom.m instead! This m file is provided to show 
% that the "standard Nystrom method" is not optimal!
% 
% Note: Use FindRep.m to find the landmark set Z. 
%
% Farhad Pourkamali-Anaraki, E-mail: Farhad.Pourkamali@colorado.edu
% University of Colorado Boulder
%
%{
Inputs:
    - X: input data matrix of size pxn, where p is the dimension and n is
    the number of samples
    - Z: landmark matrix of size pxm, where m is the number of landmark
    points
    - desiredRank: target rank
    - kernel.type: kernel type
        1) 'RBF': Gaussian - k(x,y) = exp(-gamma.*|x-y|_2^2) [parameter: gamma]
        2) 'Poly': Polynomial - k(x,y) = (x'y+c).^d [parameters: c & d]
    - kernel.par: parameters for kernels, gamma, c and d. For polynomial
    kernels, the order should be [degree d,constant c]. 

Output:
    - L: matrix of size (nxdesiredRank) which gives approximation of the
    kernel matrix K, in the form of L*L' 
%}
m = size(Z,2);
if size(X,1)~=size(Z,1), error('The given landmark set is not valid!'); end
if m < desiredRank, error('Select more landmark points!'); end


% start switch
switch kernel.type
    case 'RBF'
        gamma = kernel.par;
        W = exp( -gamma.*sqdist(Z) );   % W: m*m
        C = exp( -gamma.*sqdist(X,Z) ); % C: n*m 
    case 'Poly'
        d = kernel.par(1); c = kernel.par(2); 
        W = (Z' * Z + c).^d; % W: m*m
        C = (X' * Z + c).^d; % C: n*m 
end % end switch 

% make sure W is symmetric
W = (W + W')/2; 

% Eigenvalue Decomposition
[UW,SW] = eig(W);
[SW,I] = sort(diag(SW),'descend');
UW = UW(:, I);

SW = 1 ./ sqrt(SW(1:desiredRank));
UW = bsxfun(@times , UW(:,1:desiredRank) , SW');
L  = C * UW; % approximated by L * L'
end

