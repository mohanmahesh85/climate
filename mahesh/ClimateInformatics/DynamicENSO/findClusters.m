function [C,V] = findClusters( p2, nlat, nlon, k, ...
                            useLSMask, lsmask )

tic();
% construct affinity matrix
A = AffinityMatrix( p2, nlat, nlon, useLSMask, lsmask );
toc()

tic;
% Construct the graph Laplacian
temp = sum(A);
temp = 1./ sqrt(temp);
ii = isfinite( temp );
indexes = ii;
temp = temp(ii);

%D = diag( temp );
%L = eye(size(D)) - D*A(ii,ii)*D;
n = length(temp);
D = sparse( 1:n, 1:n, temp );
L = speye( n ) - D * A(ii,ii) * D;

toc()
tic()
% compute the eigen values/vectors of the Laplacian
    del = speye(size(L)) * 0.0001; 
    opts.issym = 1;
    [ eig_vec, eig_val ] = eigs(L + del, k,'SM', opts);
    toc
    
    tic
% determine the projection of the points on the required subspace
    %  eg_val = zeros(size(eig_val,1),1);
    %     for i = 1:size(eg_val,1)
    %         eg_val( i, 1 ) = norm(eig_val(:,i));
    %     end    
    eg_val = diag( eig_val );
    [eg_val, indices] = sortrows(eg_val);
    
    indices = indices( eg_val ~= 0 );
    
    X = zeros( size(eig_vec,1), k );
    
    for i = 1:min( k, length(indices) )
        X(:,i) = eig_vec(:,indices(i));
        
    end
    
    X = orth(X);
    
    % normalize X
    Y = X;
%     Y = zeros( size(X) );
%     for i=1:min( k, length(indices) )
%         Y(:,i) = X(:,i) / norm(X(:,i));
%     end
    

% compute the kmeans clustering for the resulting set of points
% and return the cluster assignment
[C, V] = vl_kmeans( Y', k, 'algorithm', 'elkan', 'initialization', 'plusplus') ;
toc

% reformat the variables so that the input sizes match
V2 = zeros( numel(p2), 1 );
ii = find( indexes ~= 0 );
V2(ii) = V;

%V = reshape( V2, size(p2,1), size(p2,2) );
m = size(p2,1);
n = size(p2,2);
V = zeros( m, n );

for jj = 1:numel(p2)
    j = mod( jj-1, n ) + 1;
    i = ceil( jj / n);
    V(i,j) = V2(jj);
end


end

function A = AffinityMatrix( p, lat, lon, useLSMask, lsmask )

m = size(p,1);
n = size(p,2);
A = zeros( m * n );
w = 2;

for ii = 1:m
    for jj = 1:n
        if isLand( useLSMask, lsmask, ii, jj )
            continue;
        end
        % for pixel ii,jj check its immediate neighborhood
        
        for kk = -w:w
            for ll = -w:w
                if ii + kk <= 0 || ii + kk > m ...
                        || jj + ll <= 0 || jj + ll > n...
                        || kk == 0 || ll == 0
                    continue;
                end
                
                if isLand( useLSMask, lsmask, ii+kk, jj+ll )
                    continue;
                end
                
                % use the unilateral filter
                %A( n*(ii-1) + (jj), n*(ii-1+kk) + (jj+ll) ) = ...
                %    affinity( p(ii,jj), p(ii+kk, jj+ll) ); 
                
                % use the bilateral filter
                A( n*(ii-1) + (jj), n*(ii-1+kk) + (jj+ll) ) = ...
                    affinity_b( p(ii,jj), p(ii+kk, jj+ll), ...
                    ii, jj, ii+kk, jj+ll ); 
            end
        end
    end
end

A = sparse(A);
    
end


function a = affinity( v1, v2 )
    s = 0.25;
    a = 0;
    if isnan( v1 ) || isnan( v2 )
        return;
    end
    a = exp( -0.5 * (v1 - v2)^2 / s^2 );
end

function a = affinity_b( v1, v2, x1, y1, x2, y2 )
    if v1*v1 < 0.1
        v1 = 0;
    elseif v1*v1 > 1
        v1 = sign(v1) * 1;
    end
    if v2*v2 < 0.1
        v2 = 0;
    elseif v2*v2 > 1
        v2 = sign(v2) * 1;
    end
    n = max( [0.01, v1, v2] ); 
    s = 0.15;
    s_b = 3;
    a = 0;
    if isnan( v1 ) || isnan( v2 )
        return;
    end
    a = exp( -0.5 * (v1 - v2)^2 / s^2 );% * ...
        %exp( -0.5 * ((x1 - x2)^2 + (y1-y2)^2) / s_b^2 );
end

