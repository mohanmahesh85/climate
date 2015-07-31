function [x,y] = createTexasMonthlyProblem( p, nlat, nlon, mon, t, lag  )
    
% extract the variables to be predicted

    y = extractPredictand( p, nlat, nlon );
    y = tsmovavg( y, 's', 3, 1 );
    
    x = reshape( p, size(p,1), [] );
    x = tsmovavg( x, 's', 3, 2 );
    x = reshape( x, size(p,1), size(p,2), size(p,3) );
    
    t2 = double( squeeze(t(:)) );
    t3 = datenum('1-Jan-1800');
    mm = month( t3 + t2 );
    
    ii = find( mm == mon );
    y = y(ii,:,:);
    
    mon2 = mod( mon - lag, 12 ) + 1;
    ii = find( mm == mon2 );
    
    if mon2 > mon 
        x = x(ii(1:end-1),:,:);
        y = y(2:end);
    else
        x = x(ii,:,:);
    end
    
end