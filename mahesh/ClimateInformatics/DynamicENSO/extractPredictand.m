function y = extractPredictand( x, nlat, nlon )

    t = size(x,1);
    y = zeros( t, 1 );
%     count = 0;
%     
%     w = zeros( size(x,2), size(x,3) );
%     % box centered over Texas (254-266 longitude, 26-36N latitude)
%     for ii = 1:size(x,2)
%         for jj = 1:size(x,3)
%             
%             if nlat(ii) < 35 && nlat(ii) > 27 && ...
%                     nlon(jj) <= 266 && nlon(jj) >= 252
%                 
%                 z = x(:,ii,jj);
%                 z(isnan(z)) = 0;
%                 y = y + z;
%                 w(ii,jj) = 1;
%                 count = count + 1;
%             end
%             
%         end
%     end
%     y( isnan(y) ) = 0;
%     y = y / count;

    idx1 = find( nlat < 35 & nlat > 27 );
    idx2 = find( nlon < 266 & nlon >= 252 );
    
    l1 = length(idx1);
    l2 = length(idx2);
    
    temp = zeros( t, l1 * l2 );
    for ii = 1:l1
        for jj = 1:l2
            temp(:, (ii-1)*l2 + jj) = x(:,idx1(ii), idx2(jj));
        end
    end
    
    y = nanmean( temp' )';
end
