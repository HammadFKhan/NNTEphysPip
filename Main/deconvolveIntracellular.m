%%
y = typeA;
c=exp(-(1:length(y))./80);
yc=typeB;
ydc=deconv(yc,c).*sum(c);    

figure,
subplot(2,2,1); plot(y); title('original y');
subplot(2,2,2); plot(c);title('c'); subplot(2,2,3);
plot(yc(1:length(y))); title('yc'); subplot(2,2,4);
plot(ydc);title('recovered y')
%%
figure,hold on
plot(yc,'k'),plot(ydc,'r'),ylim([-70 100])