
function dispClusters( p, nlat, nlon, v, useLSMask, lsmask )

subplot( 1, 2, 2 );
u = unique( v );
k = size( u,1 ) - 1;



% determine the extent of the figure
latlim = double( [min(nlat(:)) max(nlat(:))] );
lonlim = double( [min(nlon(:)) max(nlon(:))] );


h = worldmap( latlim, lonlim );
%set( h, 'FontSize', 0.1 );
hold on;
load coast;
plotm( lat, long, 'k' );


% here we use the average intensity over the cluster to display each
% cluster
d = zeros( size(v) );
avg_intensity = zeros( k, 1 );

for kk = 1:k
    ii = find( v == kk & ~isnan(p) );
    if isempty(ii)
        continue;
    end
    mx = mean( p(ii) );
    avg_intensity(kk) = mx(1);
    d(ii) = avg_intensity(kk);
    b = zeros( size(d) );
    b(ii) = 1;
    [B_kk, L N A] = bwboundaries( b, 'noholes' );
    if N > 1
        B(kk).b = B_kk{1};
    else
        B(kk).b = B_kk;
    end
    
%     if kk > 10 
%         figure; imagesc(b); hold on;
%         boundary = uint16( cell2mat( B_kk) );
%         plot( boundary(:,2), boundary(:,1), 'r');
%     end
end


[lat, lon] = meshgrat( nlat(1:end), nlon(1:end) );
%pcolorm( lat, lon, lsmask );
%caxis([0,k]);
%h = pcolorm( lat, lon, v );
caxis([-4,4]);
h = pcolorm( lat, lon, d );
colorbar

if useLSMask == 1
    mask = ones(size(d) );
    ii = find( v == 0 );
    mask( ii ) = 0;
    
    %set(h,'alphadata',mask,'facealpha','flat','edgecolor','none');
end

% hold on
% % draw the outlines for each of the clusters
% for kk = 1:length(B)
%     if ~isempty( B(kk).b )
%         if strcmp( class( B(kk).b ), 'double' ) == 1
%             boundary = B(kk).b;
%         else
%             boundary = uint16(cell2mat( B(kk).b ));
%         end
%         plotm(nlat(boundary(:,1)), nlon(boundary(:,2)), 'w')
%     end
% end

%linem( [-170, -120], [-5,-5] ,'k' );
% linem( [-120, -120], [-5,5] ,'k' );
% linem( [-120, -170], [5,5] ,'k' );

plotm( [-5, -5, 5, 5, -5], [-170, -120, -120,-170, -170], ...
    'magenta');

title(['Segmentation based on Temp Anomalies for k = ' num2str(k)]);

end
