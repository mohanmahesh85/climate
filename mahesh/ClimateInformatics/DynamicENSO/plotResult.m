function plotResult( folder, cm, rm, numClusters )

    path = '../data/mlost/air_mon_anom.nc';

    % read in the data
    d = ncdataset( path );
    
    % get the individual variables from the dataset
    nlat=ncvariable(d,'lat');
    nlon=ncvariable(d,'lon');
    nlat = double( nlat(:) );
    nlon = double( nlon(:) );
    fname = ['Grid_' cm '_' rm '_' num2str(numClusters) '.txt' ];
    
    grid = load( [folder fname] );
    
    plotGridMap( grid, nlat, nlon );
    title( [cm ' ' rm ' ' num2str(numClusters)] );
    colorbar
    
end