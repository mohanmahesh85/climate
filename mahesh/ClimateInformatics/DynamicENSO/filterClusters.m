function [v] = filterClusters( p, v, k )

avg_intensity = zeros( k, 1 );
for kk = 1:k
    ii = find( v == kk & ~isnan(p) );
    avg_intensity(kk) = mean( p(ii) );
end

for kk = 1:k
    ii = find( v == kk & ~isnan(p) );
    b = zeros( size(v) );
    b(ii) = 1;
    [B_kk, L N A] = bwboundaries( b, 'noholes' );
    if N > 1
        for jj = 1:N
            % check to see if clusters are too small
            if size( B_kk{jj}, 1 ) < 5
                % if yes merge them with most similar cluster
                c = mostSimilarCluster( B_kk{jj}, v, ...
                            kk, avg_intensity(kk), ...
                            avg_intensity );
               
                v = update( v, B_kk{jj}, c );
            end
            
        end
        
    end
end

end

function c = mostSimilarCluster( B, v, kk, val, arr )

    neighbours = [];
    % find the neighbouring clusters
    for ii = 1:size(B,1)
        for jj = -1:1
            for ll = -1:1
                x = B(ii,1) + jj;
                y = B(ii,2) + ll;
                if  ~(x <= 0 || x > size(v,1) || ...
                        y <= 0 || y > size(v,2) ||...
                        v(x,y) == kk || v(x,y) == 0 )
                    neighbours = [neighbours v(x, y)];
                end
            end
        end
    end
    
    neighbours = unique( neighbours );
    
    % for all the neighbours, check for the closest in intensity
    diff = (arr( neighbours ) - val).^2;
    i = find( diff == min( diff ) );
    if isempty( i ) 
        c = 0;
    else
        c = neighbours(i);
    end
end


function v = update( v, B, c )
    
    for ii = 1:size(B,1)
        v( B(ii,1), B(ii,2) ) = c;
    end
    
end