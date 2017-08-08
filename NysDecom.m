function [U , D] = NysDecom(X , Z , desiredRank , kernel)
% option 1:  [U , D] = NysDecom(X , Z , desiredRank , kernel) 
% option 2:  L = NysDecom(X , Z , desiredRank , kernel) 
%
% option 1: 
% [U , D] = NysDecom( X , Z , desiredRank , kernel ) computes eigenvalues
% and eigenvectors of the kernel matrix K.
%
% option 2: 
% L = NysDecom(X , Z , desiredRank , kernel) 
% L is a matrix of size (nxdesiredRank) which gives an approximation of the kernel matrix K,
% in the form of L * L'. 
%
% Note: Use 'FindRep.m' to find the landmark set Z.  
%
% Farhad Pourkamali-Anaraki, E-mail: Farhad.Pourkamali@colorado.edu
% University of Colorado Boulder
%{
Inputs:
    - X: input data matrix of size pxn, where p is the dimension and n is
    the number of samples
    - Z: landmark matrix of size pxm, where m is the number of landmark
    points (Use FindRep.m to find the landmark points)
    - desiredRank: target rank
    - kernel.type: kernel type
        1) 'RBF': Gaussian - k(x,y) = exp(-gamma.*|x-y|_2^2) [parameter: gamma]
        2) 'Poly': Polynomial - k(x,y) = (x'y+c).^d [parameters: c & d]
    - kernel.par: parameters for kernels, gamma, c and d. For polynomial
    kernels, the order should be [degree d,constant c]. 
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

W = (W + W')/2; % make sure W is symmetric
[Q , R] = qr (C, 0); % thin QR decomposition of matrix C
M = (R * pinv(W) * R'); M = (M + M')/2;
[V,D] = eig(M);
[D,I] = sort(diag(D),'descend');
V = V(:, I);
E = Q * V(: , 1:desiredRank);

% if eigenvalue decomposition is desired (two outputs)
if nargout > 1
    U = E; 
    D = diag(D(1:desiredRank));
else     % if low-rank approximation is desired (one output)
    U = bsxfun(@times, E, sqrt(D(1:desiredRank))');
end

end

