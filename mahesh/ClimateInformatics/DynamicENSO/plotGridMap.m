function plotGridMap( g, nlat, nlon )

% determine the extent of the figure
    latlim = double( [min(nlat(:)) max(nlat(:))] );
    lonlim = double( [min(nlon(:)) max(nlon(:))] );

    % set up the figure
    % figure('Color','white');

    
    worldmap( latlim, lonlim );
    hold on;
    load coast;

    plotm( lat, long, 'k' );
    [lat, lon] = meshgrat( nlat(1:end), nlon(1:end) );
    %caxis([0,40]);
    h = pcolorm( lat, lon, g, 'FaceAlpha', 0.7 );
    
    % draw the enso box
    plotm( [-5, -5, 5, 5, -5], [-170, -120, -120,-170, -170], ...
    'magenta');
end