% Plot time-domain causal graph, p-values and significance.
figure(); clf;
sgtitlex('Pairwise-conditional Granger causality - time domain');
% subplot(1,3,1);
plot_pw(GCoutput.F);
%%
figure,
title('Pairwise-conditional GC');
plot_pw(pval);
title(['p-values (' tstat '-test)']);
subplot(1,3,3);
plot_pw(sig);
title(['Significant at \alpha = ' num2str(alpha)]);
%% Plot spectral causal graph.
figure();
sgtitlex('Pairwise-conditional Granger causality - frequency domain');
plot_spw(GCoutput.P([5:6,46:47],[5:6,46:47],:),GCoutput.ops.fs,[4 500]);
%% find significant GC
[sigx,sigy] = find(GCoutput.sig==1);
M2M2 = [];
M1M1 = [];
M1M2 = [];
M2M1 = [];
for n = 1:length(sigy)
    if sigy(n)<=30 & sigx(n)<=30
        M2M2(n,1) = GCoutput.F(sigx(n),sigy(n));
    elseif sigy(n)<=30 & sigx(n)>30
        M1M2(n,1) = GCoutput.F(sigx(n),sigy(n));
    elseif sigy(n)>30 & sigx(n)<=30
        M2M1(n,1) = GCoutput.F(sigx(n),sigy(n));
    elseif sigy(n)>30 & sigx(n)>30
        M1M1(n,1) = GCoutput.F(sigx(n),sigy(n));
    end
end
M2M2(M2M2==0) = [];
M2M1(M2M1==0) = [];
M1M2(M1M2==0) = [];
M1M1(M1M1==0) = [];
%%
bins = 10.^(-4:0.1:.1);
figure,
histogram(M1M1,bins,'Normalization','probability','displaystyle','stairs'),set(gca, "XScale", "log"),hold on,set(gca,'TickDir','out','fontsize',16),box off
histogram(M2M2,bins,'Normalization','probability','displaystyle','stairs')
figure,
histogram(M2M1,bins,'Normalization','probability','displaystyle','stairs'),set(gca, "XScale", "log"),hold on,set(gca,'TickDir','out','fontsize',16),box off
histogram(M1M2,bins,'Normalization','probability','displaystyle','stairs')
legend('M2M1','M1M2')

