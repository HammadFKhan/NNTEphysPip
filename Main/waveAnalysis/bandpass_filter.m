function xo = bandpass_filter(x,f1,f2,Fs)
ct = [f1 f2];
ct = ct / (Fs/2);
filter_order = 4;
[b,a] = butter( filter_order, ct );
% proceed with filtering
for rr = 1:size(x,1) % electrode
    disp(['Electrode: ' num2str(rr)])
        if ~all( isnan(squeeze(x(rr,:))) )
            xo(rr,:) = filtfilt( b, a, squeeze(x(rr,:)) );
        end
end
end