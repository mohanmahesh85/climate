function [v] = clusterTimeSeries( p, k,...
                    useLSMask, lsmask, nlat, nlon, cm )
    if cm == 1
        
        if ndims(p) == 2
            p = reshape( p, size(p,1), size(lsmask,1), size(lsmask,2) );
        end
        
        if exist('A.txt','file')
            data_dump = load( 'A.txt');
            A = sparse(data_dump(:,1), data_dump(:,2), data_dump(:,3));
        else
            A = AffMat_ts( p, useLSMask, lsmask );
            [i, j, l] = find(A);
            data_dump = [i j l];
            save('A.txt','data_dump', '-ascii');
        end
        [~,v] = spectralClustering( A, k, size(p,2), size(p,3) );
        v = reshape( v, size(lsmask,1) * size(lsmask,2), 1 );
        
    elseif cm == 2
        
        v = kmeansClustering( p, k, size(p,2), size(p,3),...
           useLSMask, lsmask );
       
    end
    %v = hierarchicalClustering( p, k, size(p,2), size(p,3),...
    %            useLSMask, lsmask, nlat, nlon ); 
    
%     % try to determine points that are constantly switching between
%     % clusters and discard them
    
%     % step 1 : determine the jumpy points
%     prevV = [];
%     numTries = 1;
%     W = zeros( size(p,2) * size(p,3) );
%     for ii = 1:numTries
%         % determine the clusters
%         v = kmeansClustering( p, k, size(p,2), size(p,3),...
%             useLSMask, lsmask );
%         
%         % find the disagreements
%         if ~isempty( prevV )
%             [w,d] = disagreements( v, prevV );
%             W = W + w;
%         end
%         
%         prevV = v;
%     end
    
%     W = reshape( sum(W), size(v) );
%     Q = W;
%     Q(Q>100) = -100;
%     figure; 
%     plotGridMap( Q, nlat, nlon );
%     
%     % step 2 : discard them
%     p( :, W > 100 ) = -100;
    
    % step 3 : recluster using the refined points
    
        
end

function v = hierarchicalClustering( x, k, s1, s2,...
                useLSMask, lsmask, nlat, nlon )
    
    x(isnan(x)) = 0;
    x2 = reshape( x, size(x,1),[] );
    if useLSMask == 1
        x2(:, lsmask < 0.5 ) = 0;
    end
    
    A = pdist( x2' );
    Z = linkage(A); % load('Z.txt');%
    v = cluster( Z, 'maxclust', k );
    
    % reformat the variables so that the input sizes match
    v = reshape( v, s1, s2 );



end

function v = kmeansClustering( x, k, s1, s2, ...
                useLSMask, lsmask )

    
    x(isnan(x)) = 0;
    if ndims(x) > 2 
        x2 = reshape( x, size(x,1),[] );
    else
        x2 = x;
    end
    
    if useLSMask == 1
        x2(:, lsmask < 0.5 ) = NaN;
        %zz = std(x2);
        %x2(:, zz < 0.05 ) = NaN;
    end
    
    %[v2, c] = kmeans( x2', k, 'EmptyAction', 'singleton', 'Start', 'plus', 'Distance', 'correlation' );
    [v2, c] = kmeans( x2', k, 'EmptyAction', 'singleton', 'Start', 'plus' );
    
    v = reshape( v2, s1, s2 );
    v( isnan(v) ) = 0;
end

% assuming v1 and v2 are of the same size
function [w,d] = disagreements( v1, v2 )

    d = 0;
    v1 = reshape( v1, [], 1);
    v2 = reshape( v2, [], 1);
    
    w = zeros( size(v1) );
    
    for ii = 1:size(v1,1)
        for jj = ii+1:length(v1)
            x = xor(isSameCluster(v1,ii,jj), isSameCluster(v2,ii,jj));
            d = d + x;
            w( ii, jj ) = x;
            w( jj, ii ) = x; % maintaining a symm matrix
        end
    end
    
end

function d = isSameCluster( v, ii, jj )
    d = (v(ii) == v(jj));
end


function A = AffMat_ts( p, useLSMask, lsmask )

t = size(p,1);
m = size(p,2);
n = size(p,3);
A = zeros( m * n ) + NaN;
w = 10;

p(isnan(p)) = 0;
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
                
                
                % use the bilateral filter
                pp = [1:t];
                a = aff_ts( p(pp,ii,jj), p(pp, ii+kk, jj+ll) );
                a = exp( -0.5 * (a - 1 )^2 / 14^2 );
                d = exp( -0.5 * ( kk^2 + ll^2 ) / w );
                %A( n*(ii-1) + (jj), n*(ii-1+kk) + (jj+ll) ) = ...
                A( (ii) + m*(jj-1), (ii+kk) + m*(jj-1+ll) ) = ...
                     a * d;
            end
        end
    end
end

% use this if you are using l2
%nansum(nansum(A)) / numel(A)

%s = max( nanstd( A ) )
%A = sparse(exp(-0.5 * A.^2 / (s^2) ));
A(isnan(A)) = 0;

% use this if you are using cross corr
%A = A / max(max(A));


end

function [a,x] = aff_ts( p1, p2 )
% using a max lag of 5 for cross correlation
    a = NaN; x = 0;
%     ll1 = length(p1); ll2 = length(p2);
%     % return if you have too many nans
%     if sum( isnan( p1 ) ) > 0.9 * ll1 || ...
%             sum( isnan(p2) ) > 0.9 * ll2
%         return;
%     end
     
%     % compute the cross correlation of the time series
%     p1( isnan(p1) ) = 0;
%     p2( isnan(p2) ) = 0;
%     [r, lags] = xcorr( p1, p2, 5 );
%         a = max(r);
%         ii = find( r == a );
%         x = lags(ii(1));
     %end

     % use for l2
    p1( isnan(p1) ) = 0;%nanmean(p1);
    p2( isnan(p2) ) = 0;%nanmean(p2);
    
    if norm(p1) ~= 0 && norm(p2) ~= 0
        a = norm(p1-p2);
    end
end
