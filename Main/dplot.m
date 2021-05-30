function dplot(T,Vm,lowpass_Vm,ds,ds_lowpass,...
    secondds,secondds_lowpass)
        ax1 =  subplot(3,2,1),plot(T, Vm);hold off;
        ax1.Box = 'off';
        xlabel('sec');
        ylabel('V (mV)')
        title ('Cutoff @ \infty Hz');
        
        ax2 = subplot(3,2,2),plot(T,lowpass_Vm);
        ax2.Box = 'off';
        xlabel('sec');
        ylabel('V (mV)')
        title ('Cutoff @ 500 Hz');
        xlim([0 T(end)]);
        
        % First Derivative
        ax3 = subplot(3,2,3),plot(T(2:end),ds);
        ax3.Box = 'off';
        xlabel('sec');
        ylabel('V (mV)')
        title ('V Cutoff @ \infty Hz');
        
        ax4 = subplot(3,2,4),plot(T(2:end),ds_lowpass);
        ax4.Box = 'off';
        xlabel('sec');
        ylabel('dV (mV)')
        title ('dV Cutoff @ 500 Hz')
        
        %Second Derivative
        ax5 = subplot(3,2,5),plot(T(3:end),secondds);
        ax5.Box = 'off';
        xlabel('sec');
        ylabel('d^{2}V (mV)')
        title ('d^{2}V Cutoff @ \infty Hz');
        
        ax6 = subplot(3,2,6),plot(T(3:end),secondds_lowpass);
        ax6.Box = 'off';
        xlabel('sec');
        ylabel('d^{2}V (mV)')
        title ('d^{2}V Cutoff @ 500 Hz');
        xlim([0 T(end)]);
        linkaxes([ax1 ax2],'y');
        linkaxes([ax3 ax4],'y');
        linkaxes([ax5 ax6],'y');
        linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x');
        xlim([0 T(end)]);
end