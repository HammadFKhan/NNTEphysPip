for i = 1:29
    csdTest = CSD(t'/1E6,1024,2E-5);
end
t = mean(mLFP{55}([1:43,45:64],:,:),3);
rightSide = find(s.sorted_probe_wiring(:,2)==0);
rightSide = [rightSide;find(s.sorted_probe_wiring(:,2)==36)];
[CSD] = CustomCSD(t(rightSide,:),50,1024,0,0);

CSD2 = interp2(CSD,4);
CSD3 = smoothdata(CSD2,2,'gaussian',300);
CSD4 = smoothdata(CSD3,1,'gaussian',100);

figure,
subplot(1,4,1),imagesc(CSD),colormap(jet),colorbar
subplot(1,4,2),imagesc(CSD2),colormap(jet),colorbar
subplot(1,4,3),imagesc(CSD3),colormap(jet),colorbar
subplot(1,4,4),imagesc(CSD4),colormap(jet),colorbar

figure,imagesc(CSD4),colormap(jet),colorbar

figure,stack_plot(t(rightSide,:),1,2,1024)
