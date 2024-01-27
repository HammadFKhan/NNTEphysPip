function animateWaves(trial, Waves,saveOption,waveID)
% *SPONTANEOUS WAVES DEMO*
%
% PLOT WAVE EXAMPLES     plot specific examples of spontaneous waves based
%                          on calculated values of \rho_{\phi,d}
%
% INPUT
% x - datacube (rows,cols,timepts)
% trial - trial number
% evaluation_points - evaluation time points (cf. find_evaluation_points)
% source - putative source point
%
% OUTPUT
% animated spatiotemporal plot
%

x = rad2deg(angle(Waves(trial).p));
% x = abs(Waves(trial).p);
evaluation_points = Waves(trial).evaluationPoints;
source = Waves(trial).source;
vx = Waves(trial).vx;
vy = Waves(trial).vy;

% parameters
plot_pre_time = 0; pause_length = 0.03; 

% init
M = load( 'myMap.mat' );

if saveOption == 1
    fn = 'WaveAnimation';
    writerobj = VideoWriter([fn '.avi'],"Uncompressed AVI"); % Initialize movie file
    writerobj.FrameRate = 30;
    open(writerobj);
end

ctr = 1; % wave detections counter
for jj = 1:length(evaluation_points)
    % get start and stop time, truncating if necessary
    st = Waves(trial).waveTime(jj,1) - plot_pre_time; sp = Waves(trial).waveTime(jj,2) + plot_pre_time;
    if ( st < 1 ), st = 1; end; if ( sp > size(x,3) ), sp = size(x,3); end

    % get data to plot, shuffle if option is chosen
    x_plot = x(:,:,st:sp);
    vx_plot = vx(1,st:sp);
    vy_plot = vy(1,st:sp);
    
    t1 = st:1:sp;
    t2 = st:0.25:sp;
    xplotnew = [];
    vxplotnew = [];
    vyplotnew = [];
    for iii = 1:size(x_plot,1)
        for jjj = 1:size(x_plot,2)
            xplotnew(iii,jjj,:) = interp1(t1,squeeze(x_plot(iii,jjj,:)),t2);
        end 
    end
    vxplotnew = interp1(t1,vx_plot,t2);
    vyplotnew = interp1(t1,vy_plot,t2);
    x_plot = xplotnew;
    vx_plot = vxplotnew;
    vy_plot = vyplotnew;
    
    % Paint over nans
    for n = 1:size(x_plot,3)
        x_plot(:,:,n) = inpaint_nans(x_plot(:,:,n));
    end

    % create plot
    figure(100); 
%     title( sprintf( 'trial %d, wave example %d, 0 of %d ms', trial, ctr, size(x_plot,3) ) );
%     title( sprintf( 'Detected Wave, 0 of %d ms', size(x_plot,3) ) );
    color_range = [ min(reshape(x_plot,[],1)) max(reshape(x_plot,[],1)) ];
%     color_range = [ -80 80 ];
    h = imagesc( x_plot(:,:,1) ); hold on; axis image;
    plot( source(jj,1), source(jj,2), '.', 'markersize', 35, 'color', [.7 .7 .7] );
    h2 = quiver(source(jj,1), source(jj,2),-vx_plot(1),-vy_plot(1));
    h2.Color = 'Blue';
    h2.LineWidth = 1.5;
    h2.ShowArrowHead = 'on';
    h2.MaxHeadSize = 2;
    set( gca, 'linewidth', 3, 'xtick', [], 'ytick', [], 'fontname', 'arial', 'fontsize', 16, 'ydir', 'reverse' );
    colormap(hot);
%     colormap( M.myMap ); 
    box on; xlabel( 'electrodes' ); ylabel( 'electrodes' ); caxis( color_range )

    % create colorbar
    cb = colorbar();
    %set( cb, 'location', 'southoutside' )
    set( cb, 'location', 'east' )
    set( cb, 'position', [0.7    0.325   0.06   0.4] );
    set( get(cb,'ylabel'), 'string', 'Amplitude (in \mu V)' ); set( cb, 'linewidth', 2 )
    set( get(cb,'ylabel'), 'string', 'Phase (in deg)' ); set( cb, 'linewidth', 2 )
    if (saveOption == 1)
        writeVideo(writerobj,getframe(gcf)); %grabs current fig frame
    end

    % animate plot
    for kk = 1:size(x_plot,3)
        set( h, 'cdata', x_plot(:,:,kk) );
        set(h2, 'udata',-vx_plot(1,kk),'vdata',-vy_plot(1,kk));
        set( get(gca,'title'), 'string', ...
            sprintf( 'Detected Wave, %d of %d ms', kk, size(x_plot,3) ) )
        pause(pause_length);
        if (saveOption == 1)
            writeVideo(writerobj,getframe(gcf)); %grabs current fig frame
        end
    end
    
    % increment counter
    ctr = ctr + 1;
end
if (saveOption == 1)
    close(writerobj)
    disp('Video saved to current directory')
end
