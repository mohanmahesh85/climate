function clusteringOverlay()

    % Step 0 : Set up all the files and variables
    path = '../data/mlost/air_mon_anom.nc';
    path_mask = '../data/mlost/lsmask.nc';
    %path = '/Datasets/NARR/monolevel/prate.1979.nc';

    k = 100;
    useLSMask = 1;

%     configfile
    
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
    
    V = zeros( size(x,2), size(x,3) );
    
    numEpochs = 10;
    sizeTrain = floor( length(y) / numEpochs );
    
    for ii = 1:numEpochs
        si = (ii-1) * sizeTrain + 1;
        di = ii * sizeTrain;
        
        [v] = clusterTimeSeries( x(si:di,:,:), ...
            k, useLSMask, lsmask, nlat, nlon );
        
        figure; plotGridMap( v, nlat, nlon );
    end
    
end