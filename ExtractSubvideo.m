function ExtractSubvideo(vid_name, start_frame, end_frame)
    vid = VideoReader(vid_name);
    output_name = strcat(vid_name, '_', num2str(start_frame), '-', num2str(end_frame));
    stable = VideoWriter(output_name);
    open(stable);
    numFrames = 0;
    while hasFrame(vid)
        frame = readFrame(vid);
        numFrames = numFrames + 1;
        if numFrames >= start_frame && numFrames <= end_frame
            writeVideo(stable, frame);
        end
        if numFrames == end_frame
            break
        end
    end