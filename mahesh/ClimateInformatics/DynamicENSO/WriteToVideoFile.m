function WriteToVideoFile( workingDir )
imageNames = dir(fullfile(workingDir,'','*.png'));
imageNames = {imageNames.name}';

%Create New Video with the Image Sequence

%Construct a VideoWriter object, which creates a Motion-JPEG AVI file by default.

outputVideo = VideoWriter(fullfile(workingDir,'shuttle_out.avi'));
outputVideo.FrameRate = shuttleVideo.FrameRate;
open(outputVideo)

%Loop through the image sequence, load each image, and then write it to the video.

for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir,'images',imageNames{ii}));
   writeVideo(outputVideo,img)
end

%Finalize the video file.

close(outputVideo)

end