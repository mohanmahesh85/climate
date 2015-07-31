function t = isLand( useLSMask, lsmask, x, y )

    t = 0;
    if useLSMask == 0 
        t = 0;
    else
        if lsmask(x,y) <= 0.5
            t = 1;
        end
    end
end