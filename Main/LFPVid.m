% function ensembleVid(Ensemble,AverageImage,ROIcentroid,fn)
% if nargin<4 || strcmp('',fn)
%     fn = 'ensembleVideo';
% else
%     [~,name,~] = fileparts(fn)
%     fn = name;
% end
fn = 'LFP'

set(0,'DefaultFigureWindowStyle','normal')
clear mov
writerobj = VideoWriter([fn '.avi'],'Uncompressed AVI'); % Initialize movie file
writerobj.FrameRate = 4;
open(writerobj);
[grad,~]=colorGradient([1 1 1] ,[1 0 0 ],31);

intercount = 1000;
count = 1;
for i = 1:50
    figure(count)
    axis off
    stack_plot(flip(LFP.LFP(:,((i-1)*intercount) +1:i*intercount)),1,2)
    drawnow
    mov(count) = getframe(gcf);
    writeVideo(writerobj,mov(count));
    count = count+1;
    close all
end
close(writerobj)
disp('Video saved to current directory')
