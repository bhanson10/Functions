function create_video(F,title, frame_rate)
    writerObj = VideoWriter(title, 'MPEG-4');
    writerObj.FrameRate = frame_rate;
    open(writerObj);
    for i=1:length(F)
        frame = F(i);
        writeVideo(writerObj, frame);
    end

    close(writerObj);
end
