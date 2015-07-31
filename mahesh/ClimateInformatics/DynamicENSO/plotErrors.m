function plotErrors( folder )

    err_files = dir( [folder 'Err*.txt'] );
    
    c = {'r', 'b', 'g', 'k', 'magenta', 'cyan'};
    
    % plot the NRMSE
    figure; hold on;
    
    for ii = 1: length( err_files )
        s = err_files(ii).name;
        ss = strsplit( s, {'_', '.'} );
        leg{ii} = [ss{2} '\_' ss{3}];
        
        err = load( [folder s] );
        plot( err(:,1), err(:,2), c{ii} );
    end
    t = strrep( folder, '_', ' ' );
    t = strrep( t, '/', ' ' );
    title( t );
    legend(leg);
    
    % plot the corr
    figure; hold on;
    
    for ii = 1: length( err_files )
        s = err_files(ii).name;
        ss = strsplit( s, {'_', '.'} );
        leg{ii} = [ss{2} '\_' ss{3}];
        
        err = load( [folder s] );
        plot( err(:,1), err(:,3), c{ii} );
    end
    t = strrep( folder, '_', ' ' );
    t = strrep( t, '/', ' ' );
    title( t );
    legend(leg);
    
end