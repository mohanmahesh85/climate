function testClustering()

    % Step 0 : Set up all the files and variables
    path = '../data/mlost/air_mon_anom.nc';
    path_mask = '../data/mlost/lsmask.nc';
    %path = '/Datasets/NARR/monolevel/prate.1979.nc';

    lag = 6;
    useLSMask = 1;
    enableFiltering = 0; % 0 = none, 1 = corr based, 2 = pca based
    enableInterpolation = 0;
    numClusters = [10 20 35 50];% 60 100 150];
    
    
    
    cm = 2;
    rm = 2;
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
    % use the lag to decide what needs to be predicted and what needs to be
    % used to predict. +ve lag means predicting ahead
    
    t = ncvariable(d, 'time');
  
    [x2, y2] = createTexasProblem( x, nlat, nlon, lag );
    
    
    % separate the available data into training validation and test sets
    sizeTrain = floor(0.9 * length(y2));
    size_val = floor( 0.05 * length(y2) );
    
    train_set = 1:sizeTrain;
    val_set = sizeTrain+1 : sizeTrain+size_val;
    test_set = sizeTrain+size_val+1:length(y2);
    
    train_x = x2(train_set,:);
    train_y = y2(train_set) - mean(y2(train_set) );
    
    val_x = x2(val_set,:);
    val_y = y2( val_set ) - mean(y2(val_set));
    
    test_x = x2( test_set, : );
    test_y = y2( test_set ) - mean(y2(test_set));
    
    for kk = 1 : length( numClusters )
        numIter = 2;
        
        for jj = 1:numIter
                
                % now cluster the filtered data
                ff_train_x = train_x;
                rec_train_x = train_x;
                rec_val_x = val_x;
                
                % apply the filter
                ff_train_x( :, ff == 0 ) = NaN;
                v = clusterTimeSeries( ff_train_x, numClusters(kk), ...
                    useLSMask, lsmask, nlat, nlon, cm);
                
                           
            
            %v = zeros(size(train_x));

            %zzz = reshape( v, size(x,2), size(x,3) );
            %plotGridMap( zzz, nlat, nlon );

        end
            
           
    end
    
    
end





function [y_pred, error ] = computeError_m( testx, testy, beta )
    
%testx( isnan(testx) ) = 0;
    
    testx = [ones(length(testy),1) testx ];
    
    b = testx * beta;
    c = b - repmat( testy, 1, size(b,2) );
    error = sqrt(mean(c.^2));
    
    error = 1 - (error/std(testy));
    
    y_pred = b;
end


function v2 = generateDisplayMap( v, beta )

    v2 = zeros( size(v) );
    beta = beta( :, size(beta,2) );
    
    u = unique(v);
    for ii = 2:length(u)
        idx = find( v == u(ii) );
        v2(idx) = (beta(ii-1));
    end
    
end



