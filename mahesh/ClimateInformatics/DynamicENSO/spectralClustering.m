function [C,V] = spectralClustering( A, k, s1, s2 )

    
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

    % compute the eigen values/vectors of the Laplacian
    del = speye(size(L)) * 0.0001;
    opts.issym = 1;
    [ eig_vec, eig_val ] = eigs(L + del, k,'SM', opts);
    %[ eig_vec, eig_val ] = eigs(full(L + del), k);

    
    % determine the projection of the points on the required subspace
    %  eg_val = zeros(size(eig_val,1),1);
    %     for i = 1:size(eg_val,1)
    %         eg_val( i, 1 ) = norm(eig_val(:,i));
    %     end
    eg_val = diag( eig_val );
    [eg_val, indices] = sortrows(eg_val, 1);

    %indices = indices( eg_val ~= 0 );

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
    %[C, V] = vl_kmeans( Y', k, 'algorithm', 'elkan', 'initialization', 'plusplus') ;
    [V, C] = kmeans( Y, k, 'EmptyAction', 'singleton', 'Start', 'plus') ;
    

    % reformat the variables so that the input sizes match
    numpts = s1 * s2;
    V2 = zeros( numpts, 1 );
    ii = find( indexes ~= 0 );
    V2(ii) = V;

    %V = reshape( V2, size(p2,1), size(p2,2) );
    V = zeros( s1, s2 );

    for jj = 1:numpts
        %j = mod( jj-1, s2 ) + 1;
        %i = ceil( jj / s2 );
        i = mod( jj-1, s1 ) + 1;
        j = ceil( jj / s1 );
        V(i,j) = V2(jj);
    end

end