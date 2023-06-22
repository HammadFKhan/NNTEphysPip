function [powerCWT, fwt] = calCWTSpectogram(x,t,Fs,VoicesperOctave,flimit,plot)

%Ref - https://www.cosmos.esa.int/documents/1655127/1655136/Torrence_1998_Wavelet_Guide_BAMS.pdf/001d8327-b255-3024-a2f0-ce02e33ac98f

assert( numel(x)==numel(t), 'Number of elements in x and t dont match' );

[wt,fwt] = cwt(x,'amor',Fs,'VoicesperOctave',VoicesperOctave,'FrequencyLimits',flimit); % Calculate the CWT
powerCWT = (abs(wt).^2)/abs(var(x,1)); % Estimate power from CWT 
% The power calculated is normalzied/relative to th ewhite noise power 
if (plot==1)
    yyaxis left;imagesc(t,fwt,powerCWT);colormap('jet');set(gca,'YDir','normal');title('Wavelet based Spectogram');ylabel('Frequency (Hz)');xlabel('Time (s)');
    c=colorbar;ylabel(c, 'Relative Power to white noise','FontSize',10); 
end


