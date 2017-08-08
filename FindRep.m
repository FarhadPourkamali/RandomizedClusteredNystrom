function [ Z ] = FindRep( X , m , param )
% FindRep: Find a set of m landmark points Z=[z_1,...,z_m] for a given
% dataset X=[x_1,...,x_n]

% Farhad Pourkamali-Anaraki, E-mail: Farhad.Pourkamali@colorado.edu
% University of Colorado Boulder
%{
Inputs:
    - X: input data matrix of size pxn, where p is the dimension and n is
      the number of samples
    - m: number of desired landmark points (m<n)
    - param: a structure array in MATLAB with the specified field and values:

    - param.type: we consider various strategies to find the landmark set:
        1) 'uni-sample': uniformly sample m out of n data points

        2) 'kmeans-matlab': MATLAB implementation of K-means clustering 

        3) 'feature-extract-kmeans': randomized feature extraction
            algorithm for K-means clustering [Mahoney et al.]

    - param.dim: dimension of the reduced data for randomized K-means [Mahoney et al.]
    - param.iter: maximum number of iterations for iterative algorithms
                  (default is 10)
Outputs:
    - Z: landmark matrix of size pxm containing m representative points
%}

%{
    References: 
    [1] Zhang, Kai, Ivor W. Tsang, and James T. Kwok. "Improved Nystrom
        low-rank approximation and error analysis." In Proceedings of the
        25th international conference on Machine learning, pp. 1232-1239. ACM, 2008.

    [2] Zhang, Kai, and James T. Kwok. "Clustered Nystrom method for large scale
        manifold learning and dimension reduction." Neural Networks, IEEE Transactions
        on 21, no. 10 (2010): 1576-1587.

    [3] Arthur, David, and Sergei Vassilvitskii. "k-means++: The advantages of careful
        seeding." In Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete 
        algorithms, pp. 1027-1035. Society for Industrial and Applied Mathematics, 2007.

    [4] Boutsidis, Christos, Anastasios Zouzias, Michael W. Mahoney, and Petros Drineas. 
        "Randomized dimensionality reduction for kmeans clustering." Information Theory, 
        IEEE Transactions on 61, no. 2 (2015): 1045-1062.

    [5] Pourkamali-Anaraki, Farhad, and Stephen Becker. "Preconditioned Data Sparsification
        for Big Data with Applications to PCA and K-means." 
        arXiv preprint arXiv:1511.00152 (2015).
%}


% dataset 
[p,n] = size(X); % n is the number of data points in the dataset

if n<m, error('number of samples n must be >= number of landmark points m'); end


if ~isfield(param, 'iter')
     param.iter = 10; % default value 
end

% start switch
switch param.type
 
    case 'uni-sample'
         indx = randsample(n,m);
         Z    = X(:,indx);
         
 
    case 'kmeans-matlab'
        [~,Z] = kmeans(X',m,'Replicates',1,'MaxIter',param.iter); % m centroids
        Z     = Z';
        
    case 'feature-extract-kmeans' 
        r     = param.dim;
        R     = randsrc(r,p,[-1 1]);
        Y     = R * X;
        [idx,~] = kmeans(Y',m,'Replicates',1,'MaxIter',param.iter); % m centroids
        Z = computeMeans( X, m, int32(idx'-1) ); % one more pass over the data to find cluster centers
      
    otherwise % otherwise
        error('Please change param.type to uni-sample or kmeans-matlab or feature-extract-kmeans!')
end % end switch

end

