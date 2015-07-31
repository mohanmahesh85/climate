function [x,y] = createTexasProblem( p, nlat, nlon, lag  )
    
% extract the variables to be predicted

    y = extractPredictand( p, nlat, nlon );
    y = tsmovavg( y, 's', 3, 1 );
    
    x = reshape( p, size(p,1), [] );
%     x = tsmovavg( x, 's', 3, 2 );
    x = reshape( x, size(p,1), size(p,2), size(p,3) );
    
    
    y = y(lag+1:end);
    x = x(1:end-lag,:,:);

end