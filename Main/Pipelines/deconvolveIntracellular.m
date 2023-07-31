%%
y = typeA(1:10000);
yc = typeB;
c=exp(-(1:length(y))./80);
ydc=deconv(yc,c).*sum(c);    

figure,
subplot(2,2,1); plot(y); title('original y');
subplot(2,2,2); plot(c);title('c'); subplot(2,2,3);
plot(yc); title('yc'); subplot(2,2,4);
plot(ydc);title('recovered y')
%%
Fc = [280];
Wn = Fc./(Fs/2);
b = butter(1,Wn,'high');
dat = filter(b,1,typeB);
%%
figure,hold on
plot(yc,'k'),plot(ydc,'r'),ylim([-70 100])
%%
data = filter(cell2mat(Num),cell2mat(Den),typeB);
xdata = 0:1/Fs:(length(data)-1)/Fs;
figure,plot(0:1/Fs:(length(data)-1)/Fs,data*80000)
hold on
plot(0:1/Fs:(length(data)-1)/Fs,typeB,'r')
