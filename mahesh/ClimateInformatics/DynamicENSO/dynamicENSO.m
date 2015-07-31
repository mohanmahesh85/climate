function dynamicENSO()
path = '../data/mlost/air_mon_anom.nc';
path_mask = '../data/mlost/lsmask.nc';
%path = '/Datasets/NARR/monolevel/prate.1979.nc';

k = 30;
useLSMask = 1;

% read in the data
d = ncdataset( path );
lsmask = ncdataset( path_mask );


% get the individual variables from the dataset
p = ncvariable( d, 'air');
nlat=ncvariable(d,'lat');
nlon=ncvariable(d,'lon');
nlat = double( nlat(:) );
nlon = double( nlon(:) );
ss = p.size;

%month = repmat( [1:12]', [ss(1)/12 1] );
% year = 1880 + [0:ss(1)/12];

t = ncvariable(d, 'time');
t2 = double( squeeze(t(:)) );
t3 = datenum('1-Jan-1800');

ll = ncvariable( lsmask, 'lsmask');
lsmask = double( squeeze( ll(:,:,:) ) );

% determine the extent of the figure
latlim = double( [min(nlat(:)) max(nlat(:))] );
lonlim = double( [min(nlon(:)) max(nlon(:))] );

% set up the figure
figure('Color','white');

global mei_data;
m = mei(12, 1970);

for ii = 1:ss(1)
    subplot( 1, 2, 1);
    worldmap( latlim, lonlim );
    hold on;
    load coast;
    
    plotm( lat, long, 'k' );
  
   
    %year = 1880 + ceil( double(ii) / 12 );
    mm = month( t3 + t2(ii) );
    yy = year( t3 + t2(ii) );
    
    title( ['SST Anomaly for Month = ' num2str(mm) ' Year = ' ...
        num2str(yy)] );
        
    % prepare the grid variables for display
    [lat, lon] = meshgrat( nlat(1:end), nlon(1:end) );
    p2 = double( squeeze(p(ii,1:end,1:end) ));
    
    %display the deviations
    caxis( [-4,4] );
    pcolorm(lat, lon, p2 );
    colorbar;
   
    
    % plot the value of current mei
%     m = [];
%     for jj = 5:-1:1
%         mm = month( t3 + t2(ii -jj) );
%         yy = year( t3 + t2(ii - jj) );
%         m = [m mei(mm,yy)];
%     end
%     
%     subplot( 1,2,2 );
%     plot( 1:length(m), m );
    
    
    
    % find clusters here
    % p = ncvariable( d, 'air');

    % convert these into double arrays and vectors    
    p2 = double( squeeze( p(ii,:,:) ) );

    [c, v] = findClusters( p2, nlat, nlon, k, ...
                useLSMask, lsmask );
    %v = vl_slic( single(p2), 10, 0.01 );
    
    % Just merge clusters which are too small
    %[v] = filterClusters( p2, v, k );
    
    % display the clusters here
    dispClusters( p2, nlat, nlon, v, useLSMask, lsmask );
    
    hold off;
    tightmap
    pause(0.01);
    
    print( ['spec_nodist/' num2str(yy) '_' num2str(mm, '%0.2d') ], '-dpng' ); 
    %print( ['' num2str(yy) '_' num2str(mm, '%0.2d') ], '-dpng' ); 
    %waitforbuttonpress
end
end

