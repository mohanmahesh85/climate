% This function takes two grids arr1, arr2 of different resolutions
% and and scales arr2 to return arr3 which is at the same resolution
% as arr1
% use this only when arr2 is at a finer resolution than arr1
function arr3 = make_same_size( arr1, lat1, lon1,...
                         arr2, lat2, lon2 )

% arr3 = zeros( size(arr1) );
% 
% for ii = 1:length(lat1)
%     for jj = 1:length(lon1)
%         arr3(ii,jj) = comp_av_val( lat1(ii), lon1(jj),...
%                                     arr2, lat2, lon2 );
%     end
% end

arr3 = imresize( arr1, length(lat2), length(lon2) );

end

% compute the average value for the grid cell at lat, lon
% based on the values in arr
function av = comp_av_val( lat, lon, ...
                            arr, lat2, lon2 )
                        
    
end