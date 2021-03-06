% This function reads in the  Multivariate ENSO index from 
% mei.txt and mei.ext.txt. 
% Ref : Wolter, K., and M. S. Timlin, 2011: El Niño/Southern Oscillation
% behaviour since 1871 as diagnosed in an extended multivariate ENSO index (MEI.ext). Intl. J. Climatology, 31
% link : http://www.esrl.noaa.gov/psd/enso/mei.ext/#data
function m = mei( month, year )

global mei_data;

if( size( mei_data, 1 ) == 0 )
    m1 = load( '../data/mei.txt' );
    m2 = load( '../data/mei.ext.txt' );
    
    % combine the two. While combining pref given to m1
    fy  = m1(1,1);
    fy2 = m2(1,1);
    
    mei_data = zeros( (fy - fy2) + size(m1,1), size(m1,2) );
    
    mei_data( 1:(fy-fy2), :)   = m2(1:(fy-fy2), :);
    mei_data( fy-fy2+1:end, :) = m1;
    
end

fy = mei_data(1,1);
if year < fy
    m = 0;
else
    m = mei_data( year - fy + 1, month + 1 );
end

end