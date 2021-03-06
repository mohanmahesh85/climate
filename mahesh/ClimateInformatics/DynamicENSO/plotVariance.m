function plotVariance()

path = '../data/mlost/air_mon_anom.nc';
path_mask = '../data/mlost/lsmask.nc';

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

ll = ncvariable( lsmask, 'lsmask');
lsmask = double( squeeze( ll(:,:,:) ) );
ii = find( lsmask > 0.5 );
lsmask =  ones( size(lsmask) );
lsmask( ii ) = 0;

latlim = double( [min(nlat(:)) max(nlat(:))] );
lonlim = double( [min(nlon(:)) max(nlon(:))] );

% plot the figure
worldmap( latlim, lonlim );
hold on;
load coast;

plotm( lat, long, 'k' );

p = double( p(:,:,:) );
[~,m,n] = size( p );

v = zeros( m, n );


for ii = 1:m
    for jj = 1:n
        if lsmask( ii,jj ) == 1
            continue;
        end
        p2 = double( squeeze(p( :, ii, jj )) );
        temp = var( p2 );
        v(ii,jj) = temp;
    end
end

ii = find( isnan(v) );
v(ii) = 0;


[lat, lon] = meshgrat( nlat(1:end), nlon(1:end) );
colormap(jet(100));
%caxis( [0, 1]);
pcolorm(lat, lon, sqrt(v) );

colorbar;

size(p)
end