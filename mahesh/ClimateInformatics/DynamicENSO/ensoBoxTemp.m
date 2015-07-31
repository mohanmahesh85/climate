function y = ensoBoxTemp( x, nlat, nlon )
%function t = ensoBoxTemp()

%plotm( [-5, -5, 5, 5, -5], [-170, -120, -120,-170, -170], ...
%    'magenta');

t = size(x,1);
    y = zeros( t, 1 );
    count = 0;
    
    w = zeros( size(x,2), size(x,3) );
    z = [];
    % box centered over Texas (254-266 longitude, 26-36N latitude)
    for ii = 1:size(x,2)
        for jj = 1:size(x,3)
            
            if nlat(ii) < 5 && nlat(ii) >= -5 && ...
                    nlon(jj) <= 266 && nlon(jj) >= 252
                
                z = [z x(:,ii,jj)];
                %z(isnan(z)) = 0;
                %y = y + z;
                w(ii,jj) = 1;
                count = count + 1;
            end
            
        end
    end
    y = nanmean( z, 2 );
    y( isnan(y) ) = 0;
end

