% This function aims to find the predictive clusters by using a
% hierarchical approach 
% 1. Group together grid cells based on the affinities between their
% associated time series
% 2. Use the aggregated grid cells to determine new features
% 3. Use new feratures for regression to predict the Texas temps
% 4. If possible use a sliding window
function dynamicENSOMonthly()

    % Step 0 : Set up all the files and variables
    path = '../data/mlost/air_mon_anom.nc';
    path_mask = '../data/mlost/lsmask.nc';
    %path = '/Datasets/NARR/monolevel/prate.1979.nc';

    
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
    % use the lag to decide what needs to be predicted and what needs to be
    % used to predict. +ve lag means predicting ahead
    lag = 3;
    t = ncvariable(d, 'time');
    
    results_folder = {'Results_consolidated2_Feb/', 'Results_consolidated2_May/', ...
        'Results_consolidated2_Aug/', 'Results_consolidated2_Nov/'};
    c_methods = {'spectral', 'kmeans'};
    r_methods = {'ridge', 'lasso', 'svr'};
    
    
    for jj = 3:length(results_folder)
        [x2, y2] = createTexasMonthlyProblem( x, nlat, nlon, (jj-1)*3+2, t, lag );

        for cm = 1:length(c_methods)
            for rm = 2:2%length(r_methods)

                numClusters = [10 20 30 40];
                numIter = 5;
                v_mean_err = zeros( length(numClusters), 1);
                v_mean_var = zeros( length(numClusters), 1);

                for kk = 1:numel(numClusters)

                    beta_mean = zeros( size(x,2), size(x,3) );
                    v_err = zeros( numIter, 1 );
                    v_var = zeros( numIter, 1 );

                    maxError = 0;
                    
                    for ii = 1:numIter
                        % first find the clusters for feature reduction
                        [v] = clusterTimeSeries( x, ...
                            numClusters(kk), useLSMask, lsmask, nlat, nlon, cm );

                        % next regress on the reduced feature set
                        [b, v_err(ii), v_var(ii)] = regressTimeSeries( x2, y2, v, rm );
                        v_err(ii)
                        % now display the regression results on the friggin map
                        %v2 = generateDisplayMap( v, b );
                        %beta_mean = beta_mean + v2;
                        
                        if v_err(ii) > maxError
                            maxError = v_err(ii);
                            % now display the regression results on the friggin map
                            v2 = generateDisplayMap( v, b );
                            beta_mean = v2;
                        end
                    end
                    %beta_mean = beta_mean / numIter;
                    v_mean_err(kk) = max( v_err );
                    idx = find( v_err == max(v_err) );
                    v_mean_var(kk) = v_var( idx(1) );
                    figure;
                    plotGridMap( beta_mean, nlat, nlon ); drawnow
                    res_file = [ results_folder{jj} 'Grid_' c_methods{cm} '_' r_methods{rm} '_' num2str(numClusters(kk)) '.txt'];
                    save( res_file, 'beta_mean', '-ascii' );
                end

                % save the errors
                res_file = [results_folder{jj} 'Err_' c_methods{cm} '_' r_methods{rm} '.txt'];
                res = [numClusters' v_mean_err v_mean_var];
                save( res_file, 'res', '-ascii' );
                v_mean_err'
            end
        end
    end
    
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


