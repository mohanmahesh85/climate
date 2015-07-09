function fig3()
  
  path = '../data/mlost/air_mon_anom.nc';
  path_precip = '../data/precip/precip.mon.mean.1x1.nc';
  
  
  % box centered over Texas (254-266 longitude, 26-36N latitude)
  latlim2 = [25 40];
  lonlim2 = [250 270];
  
  ax = worldmap(latlim2, lonlim2);
  %hold on; load coast;
  %plotm( lat, long, 'k' );

  % display the map of Texas
  states = shaperead('usastatelo', 'UseGeoCoords', true);
  faceColors = makesymbolspec('Polygon',...
    {'INDEX', [1 numel(states)], 'FaceColor', ...
    polcmap(numel(states))}); % NOTE - colors are random

  geoshow(ax, states, 'DisplayType', 'polygon', ...
      'SymbolSpec', faceColors)


  % read in the temperature anomaly data
  d = ncdataset( path );
  precip_ds = ncdataset( path_precip );
  
  
  % get the individual variables from the dataset
  p = ncvariable( d, 'air');
  nlat=ncvariable(d,'lat');
  nlon=ncvariable(d,'lon');
  nlat = double( nlat(:) );
  nlon = double( nlon(:) );
  ss = p.size;
  
  t = ncvariable(d, 'time');
  t2 = double( squeeze(t(:)) );
  t3 = datenum('1-Jan-1800');
  
  precip    = ncvariable( precip_ds, 'precip' );
  nlat2     = ncvariable( precip_ds,'lat');
  nlon2     = ncvariable( precip_ds,'lon');
  nlat2     = double( nlat2(:) );
  nlon2     = double( nlon2(:) );

  

  p_avg = zeros( ss(2), ss(3) );
  num_cases = 0;
  for ii = 1570:ss(1)

    if shouldInclude( t3 + t2(ii) ) == 0
      continue;
    end

    p2 = double( squeeze(p(ii,1:end,1:end) ));
    p3 = double( squeeze(precip(ii,1:end,1:end) ));
    precip2 = make_same_size( p2, nlat, nlon, precip, nlat2, nlon2 );
    p_avg = p_avg + p2;
    num_cases = num_cases + 1;
  end


    % prepare the grid variables for display
    [lat, lon] = meshgrat( nlat(1:end), nlon(1:end) );
    
    %display the deviations
    caxis( [-4,4] );
    pcolorm(lat, lon, p_avg/num_cases );

    drawnow;
end

% stuff is included only for June, July, August 2011
function b = shouldInclude( t )
  
  yy = year( t )
  mm = month( t );

  %if strcmp( yy, '2011' ) == 1
  %  if strcmp( mm, 'June' ) == 1 || ...
  %      strcmp( mm, 'July' ) == 1 || ...
  %       strcmp( mm, 'August' ) == 1 
  %      b = 1;
  %    end
  %  end
  
  if yy == 2011 && (mm == 6 || mm == 7 || mm == 8 )
    b = 1;
  else
    b = 0;
  end
end
