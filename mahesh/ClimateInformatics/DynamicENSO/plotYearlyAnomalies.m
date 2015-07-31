function plotYearlyAnomalies()

    % Step 0 : Set up all the files and variables
    path = '../data/mlost/air_mon_anom.nc';
    path_mask = '../data/mlost/lsmask.nc';
    %path = '/Datasets/NARR/monolevel/prate.1979.nc';
    
    % read in the data
    d = ncdataset( path );
    lsmask = ncdataset( path_mask );


    % get the individual variables from the dataset
    p = ncvariable( d, 'air');
    nlat=ncvariable(d,'lat');
    nlon=ncvariable(d,'lon');
    nlat = double( nlat(:) );
    nlon = double( nlon(:) );

    
    ll = ncvariable( lsmask, 'lsmask');
    lsmask = double( squeeze( ll(:,:,:) ) );

    
    % Step 1 : Find the aggregated grid cells
    x = double( squeeze(p(:,:,:)));
    
    % Step 2 : Set up the regression problem
    % extract the variables to be predicted
    y = extractPredictand( x, nlat, nlon );
    
    my = tsmovavg( y, 's', 3, 1 );
    figure;
    hold on;
    si = 1;
    for ii = si:12:length(my)
        dy = y(ii:ii+12-1);
        plot(1:length(dy), dy, 'r.');
        
    end

end