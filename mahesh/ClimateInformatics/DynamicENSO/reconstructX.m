function x2 = reconstructX( x, v )

    if norm(v) == 0
        x2 = x;
        return;
    end
    
    u = unique(v);
    t = size( x, 1 );

    x2 = zeros( t, length(u)-1 );

    kk = 1;
    for ii = 2:length(u)
        l = u(ii);

        % select all the grid cells belonging to this group
        [a, b] = find( v == l );


        %[kk length(a) length(u)]
        xx = zeros( t, length(a) );
        for jj = 1: length(a)
            z = x( :, a(jj), b(jj) );
            xx( : , jj) =  z;
        end
        x2(:,kk) = nanmean(xx,2);
        %x2(:,kk) = max(xx,[],2);
        kk = kk + 1;    

    end

end