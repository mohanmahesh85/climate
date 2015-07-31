function [beta, error, s] = regressTimeSeries( x, y, v, rm )

    sizeTrain = floor(0.85 * length(y));
    size_val = floor( 0.1 * length(y) );
    
    if norm(v) ~= 0 
        x2 = reconstructX( x, v );
    else
        x2 = reshape( x, size(x,1), [] );
    end
    x2(isnan(x2)) = 0;
    %x2 = t_normalize(x2);
    y2 = y;%center(y);
    train_set = 1:sizeTrain;
    val_set = sizeTrain+1 : length(y);
    test_set = sizeTrain+size_val+1:length(y);
    
    %trainx = t_normalize( x2(train_set,:) );% ones(length(train_set),1)];
    trainx = x2(train_set,:);
    trainy = y2(train_set); %t_center( y2(train_set) ); 
    
    %testx = t_normalize( x2(val_set,:));% ones(length(val_set),1)];
    testx = x2(val_set,:);
    testy = y2(val_set); %t_center( y2( val_set ) );
    
    

    if rm == 1
        beta = regRidge( trainx, trainy );
        %[train_beta, train_err, train_s] = computeError( trainx, trainy, beta );
        
        [beta, error, s] = computeError( trainx, trainy, beta );
    elseif rm == 2
        beta = regLasso( trainx, trainy );
        %[train_beta, train_err, train_s] = computeError( trainx, trainy, beta );
        
        [beta, error, s] = computeError( trainx, trainy, beta );
        %beta = beta(1:end-1,:);
    elseif rm == 3
        [trainx, trainy, testx, testy] = scaleData( trainx, trainy, ...
                                testx, testy );
        model = svmtrain( trainy, trainx, '-s 4 -t 1 -d 3' );
        y_hat = svmpredict( testy, testx, model );
        error = sqrt( mean( (testy-y_hat).^2 ) );
        error = 1 - error/std(testy);
        beta = model.SVs' * model.sv_coef;
        s = var( (testy-y_hat).^2 );
    elseif rm == 4
        % copula granger
        cause = copulaGranger( [trainx trainy]' );
        beta = cause(:,end);
        [beta, error, s] = computeError( testx, testy, beta(1:end-1) );
    end
    
end

function [trainx, trainy, testx, testy] = scaleData( trainx, trainy, ...
                                testx, testy )
    trn_data.X = trainx;
    trn_data.y = trainy;
    
    tst_data.X = testx;
    tst_data.y = testy;
    
    [trn_data, tst_data, jn2] = scaleSVM( trn_data, tst_data, trn_data, -1, 1 );
    
    trainx = trn_data.X;
    trainy = trn_data.y;
    
    testx = tst_data.X;
    testy = tst_data.y;
end
                            
function [beta] = regLasso( x, y )
    
    [beta, info ] = lasso( x, y );
    beta = [info.Intercept; beta];
    
end

function beta = regRidge( x, y )

    [beta] = ridge( y, x, 0:0.1:1, 0 );
    
end

function [beta, error, s] = computeError( testx, testy, beta )
    
    testx = [ones(length(testy),1) testx ];
    %testx( isnan(testx) ) = 0;
    
    b = testx * beta; %beta' * testx';
    %c = b - repmat( testy, 1, size(b,2) );
    c = bsxfun(@minus,testy,b);
    error = sqrt(mean(c.^2));
    s = var(c);
    
    % find the model with the lowest error
    %ii = find( error == min(error) );
    %beta = beta( :, ii(1) );
    %s = var(c(ii(1),:));
    %error = min(error);
    error = 1 - (error/std(testy));
    
    
end




function x = t_normalize( x )

    n = zeros( size( x,1 ), 1 );
    for ii = 1:size(x,1)
        n(ii) = norm( x(ii,:) );
    end
    
    for ii = 1:length(n)
        x(ii,:) = x(ii,:) / n(ii);
    end
end

function y2 = t_center( y )

    y2 = y - mean(y);
    
end
