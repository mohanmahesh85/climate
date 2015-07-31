function playAnimation( folder, cm, rm, numClusters )

    
    for kk = 1:length(numClusters)
        plotResult( folder, cm, rm, numClusters(kk) );
        waitforbuttonpress
    end
end