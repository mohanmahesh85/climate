function temporaryRun()

    % Step 0 : Set up all the files and variables
    path = '../data/mlost/air_mon_anom.nc';
    path_mask = '../data/mlost/lsmask.nc';
    %path = '/Datasets/NARR/monolevel/prate.1979.nc';

    lag = 6;
    useLSMask = 1;
    enableFiltering = 0; % 0 = none, 1 = corr based, 2 = pca based
    enableInterpolation = 0;
    numClusters = [10 20 35 50];% 60 100 150];
    
    results_folder = 'Results_consolidated_annual_6mon4/';
    c_methods = {'spectral', 'kmeans'};
    r_methods = {'ridge', 'lasso', 'svr'};
    
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
    
    %results_folder = {'Results_consolidated2_Feb/', 'Results_consolidated2_May/', ...
    %    'Results_consolidated2_Aug/', 'Results_consolidated2_Nov/'};
    
     
  
    [x2, y2] = createTexasProblem( x, nlat, nlon, lag );
    
    if enableInterpolation == 1 
        x2 = interpolateNans( x2 );
        y2 = interpolateNans( y2 );
    else
        x2( isnan(x2) ) = 0;
        y2( isnan(y2) ) = 0;
    end
    
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
    
    % first filter out the points that have a low correlation with the
    % target.
    if enableFiltering == 1
        % filter out grid cells based on correlation coeffs
        size_filter_data = floor( 0.1 * sizeTrain);
        filter_x = train_x;%( 1:size_filter_data, : );
        filter_y = train_y;%( 1:size_filter_data );
        %train_x = train_x( size_filter_data+1:end, : );
        %train_y = train_y( size_filter_data+1:end );
        
        [cc,lags] = compute_corr( filter_x, filter_y, lsmask, useLSMask );
        
        thr = 0.75;
        ff =  filter_coeff( cc, lags, lag, thr );
    elseif enableFiltering == 2
        % filter out dimensions based on the pca
        filter_x = train_x;
        coeff = pca( filter_x );
        
        ff = ( lsmask > 0.5 );
        ff = reshape( single(ff), 1, [] );
    else
        %ff = ones( size( train_x, 2 ), 1 );
        ff = ( lsmask > 0.5 );
        ff = reshape( single(ff), 1, [] );
    end
    % uncomment to visualize the filter coeff
%     zzz = reshape( ff, size(x,2), size(x,3) );
%     figure; plotGridMap( zzz, nlat, nlon );
    
      

    % finally compute the test error and other metrics for the enso temps
      x_enso = ensoBoxTemp( x, nlat, nlon );
      x_enso_test = x_enso( test_set );
      err_enso = computeTestError( x_enso_test, test_y );
      corr_enso = corr( x_enso_test, test_y );
    
    
    
    
    
    for kk = 1 : length( numClusters )
        numIter = 10;
        best_beta = zeros( numClusters(kk)+1, 1 );
        best_val_error = -100000;
        best_v = zeros( size(x,2), size(x,3) );
        mean_beta = zeros( numIter, size(best_v,1), size( best_v,2 ) );
        
        for jj = 1:numIter
            if enableFiltering == 2
                rec_train_x     = train_x * coeff(:,1:numClusters(kk));
                rec_val_x   = val_x * coeff(:,1:numClusters(kk));
                ff_train_x = rec_train_x;
                v = 0;
            else
                
                % now cluster the filtered data
                ff_train_x = train_x;
                rec_train_x = train_x;
                rec_val_x = val_x;
                
                % apply the filter
                ff_train_x( :, ff == 0 ) = NaN;
                v = clusterTimeSeries( ff_train_x, numClusters(kk), ...
                    useLSMask, lsmask, nlat, nlon, cm);
                
                
                rec_val_x   = reconstructX( rec_val_x, v );

            end
            
            
            
            %v = zeros(size(train_x));

            %zzz = reshape( v, size(x,2), size(x,3) );
            %plotGridMap( zzz, nlat, nlon );
                
            % compute the regressor for the filtered data
            [beta, train_error s] = regressTimeSeries( rec_train_x, train_y, v, rm );
            %[kk jj]
            %train_error

            
            
            % compute the training and validation error
            [y_pred, val_error] = computeError_m( rec_val_x, val_y, beta );
            idx = find( val_error == max(val_error) );
            
            
            [ val_error( idx(1) ) train_error( idx(1) ) numClusters(kk)]
            
            if max(val_error) > best_val_error
                best_beta = beta(:,idx(1));
                best_v = v;
                best_val_error = max(val_error);
                
%                 zzz = generateDisplayMap( v, best_beta );
%                 zzz = reshape( zzz, size(x,2), size(x,3) );
%                 figure; plotGridMap( zzz, nlat, nlon );
            end
            
            if enableFiltering ~= 2
                zzz = generateDisplayMap( v, beta(:,idx(1)) );
                zzz = reshape( zzz, size(x,2), size(x,3) );
                mean_beta(jj,:,:) = (zzz > 0.1 ) + -1 * (zzz < -0.1);
                bb(kk).error(jj) = max(val_error);
            end
        end
        if enableFiltering ~= 2
            bb(kk).best_beta = best_beta;
            bb(kk).best_val_error = best_val_error;
            bb(kk).best_v = best_v;
            
            % trying to average beta smartly
            s = std( bb(kk).error );
            bb(kk).mean_beta = zeros( size(mean_beta,2), size(mean_beta,3) );
            w = zeros( numIter, 1 );
            for jj = 1:numIter
                w(jj) = exp( -0.5 * (bb(kk).best_val_error - bb(kk).error(jj)) / s );
            end
            w = w / sum(w);
            for jj = 1:numIter
                bb(kk).mean_beta = bb(kk).mean_beta + w(jj) * squeeze(mean_beta(jj,:,:));
            end
            
            % display
            zzz = generateDisplayMap( bb(kk).best_v, bb(kk).best_beta );
            zzz = reshape( zzz, size(x,2), size(x,3) );
            figure; plotGridMap( zzz, nlat, nlon );
            figure; plotGridMap( bb(kk).mean_beta, nlat, nlon );
        end
    end
    
      err_method = zeros( length(numClusters), 1 );
      corr_method = zeros( length(numClusters), 1 );
      for kk = 1:length(numClusters)
         
          rec_test_x = reconstructX( test_x, bb(kk).best_v );
          [y_pred, err_method(kk)] = computeError_m( rec_test_x, test_y, bb(kk).best_beta );
          corr_method(kk) = corr( y_pred, test_y );
          
          %display
          zzz = generateDisplayMap( bb(kk).best_v, bb(kk).best_beta );
          zzz = reshape( zzz, size(x,2), size(x,3) );
          figure; plotGridMap( zzz, nlat, nlon );
          
          res_file = [ results_folder 'Grid_' c_methods{cm} '_' r_methods{rm} '_' num2str(numClusters(kk)) '.txt'];
          save( res_file, 'zzz', '-ascii' );
      end
  
    % save the errors
    res_file = [results_folder 'Err_' c_methods{cm} '_' r_methods{rm} '.txt'];
    res = [numClusters' err_method corr_method];
    save( res_file, 'res', '-ascii' );
end

function err = computeTestError( y_hat, y )

    rmse = sqrt( mean( (y_hat-y).^2 ) );
    err = 1 - rmse/std(y);
    
end


function x2 = interpolateNans( x2 )
    
    x2 = reshape( x2, size(x2,1), [] );
    
    for ii = 1:size(x2, 2)
        xx = x2(:,ii);
        
        jj = find( isnan(xx) );
        kk = find( ~isnan(xx) );
        if length(jj) > length(kk) 
            x2(:,ii) = 0;
        elseif ~isempty(jj)  
            xx(jj) = interp1( kk, xx(kk), jj, 'spline');
            x2(:,ii) = xx;
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

function ff = filter_coeff( cc, lags, lag, thr )
    
    lags( lags > 12 ) = NaN;
    lags( lags < -12 ) = NaN;
    
    sd_cc   = nanstd(abs(cc));
    sd_lags = nanstd(lags);
    
    m = nanmean( abs(cc) );
    %ff = exp( -0.5 * (abs(cc)-m/2).^2 / sd_cc^2 );% .* exp( -0.5 * (lags-lag).^2 / sd_lags^2 );
    %ff = double( abs(cc) > m ); % .* exp( -0.5 * (abs(lags)-lag).^2 / (2*sd_lags)^2 );
    %ff(isnan(ff)) = 0;
    %ff( abs(ff) < thr ) = 0;
    ff = double( abs(cc) > 50 );
end

function [c, lags] = compute_corr( tx, ty, lsmask, useLSMask )

    lags = zeros( size( tx, 2 ), size( tx, 3 ) );
    c = zeros( size( tx, 2 ), size( tx, 3 ) );
    
    
    for ii = 1:size(tx, 2)
        for jj = 1:size(tx,3)
            if useLSMask == 1 && lsmask(ii) <= 0.5 
                continue;
            end
                
            r               = tx(:, ii, jj);
            r( isnan(r) )   = 0;
            ty( isnan(ty) ) = 0;
            
            [cc, ll]        = xcorr( r, ty, 1 );
            [~, kk]          = max( abs( cc ) );
            c(ii,jj)        = cc( kk(1) );
            
            if length(kk) ~= 0 
                lags(ii,jj) = ll( kk(1) );
            end
        end
    end
    
    
end

function y2 = t_center( y )

    y2 = y - mean(y);
    
end
