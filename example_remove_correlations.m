function example_remove_correlations

%********* sample code to remove slow trends and trial averages

T = 40;  % number of trials
B = 10;  % number of time bins
Binsize = 100; % time in ms

%****** time trend ************
x = 1:B;
mrate = 1 + exp( -(0.25*x));  % exponential decay from onset response
xt = 1:T
trate = 1.5 - (xt/T);   % linear drop in rate over session
baserate = 0.05;  % prob of spiking per ms, this would be 50hz  

%***** generate data with trends ************
spcnts = cell(1,2);
for k = 1:2
   spcnts{1,k} = zeros(T,B); 
   for  ib = 1:B
       for it = 1:T
          rate = mrate(ib)*trate(it)*baserate;
          cnt = 0;
          for i = 1:Binsize
             if (rand < rate)
                 cnt = cnt + 1;
             end
          end
          spcnts{1,k}(it,ib) = cnt;
       end
   end  
end

%********** plot sample rates ************
figure;
subplot(3,3,1);
imagesc(spcnts{1,1});
title('raw cnts unit 1');
subplot(3,3,2);
imagesc(spcnts{1,2});
title('raw cnts unit 2');
subplot(3,3,3);
s1 = reshape(spcnts{1,1},[1 (T*B)]);
s2 = reshape(spcnts{1,2},[1 (T*B)]);
[r,p] = corrcoef(s1,s2);
plot(s1,s2,'ko');
title(sprintf('Raw Corr r=%+5.3f (p=%6.4f)',r(1,2),p(1,2)));

%*********** now implement normalization for slow and trial locked trends
for k = 1:2
   %********* first normalize out time locked trial trend
   bmean = mean( spcnts{1,k} );  % 1 to B values
   for i = 1:B
       spcnts{1,k}(:,i) = spcnts{1,k}(:,i) - bmean(i);
   end
   
   
   %******** second, compute Gauss weighted smoothing over trials
   tmean = mean( spcnts{1,k}' );  % 1 to T values
   gwin = 5;
   for it = 1:T
       ia = it - (2*gwin);
       if (ia < 1)
           ia = 1;
       end
       ib = it + (2*gwin);
       if (ib > T)
           ib = T;
       end
       sum = 0;
       norm = 0;
       for z = ia:ib
           weight = exp( -0.5*((it-z)*(it-z))/(gwin*gwin));
           sum = sum + (tmean(z)*weight);
           norm = norm + weight;
       end
       sum = sum / norm;
       gtmean(it) = sum;
   end
      %****** and subtract Gauss smoothed mean
   for i = 1:T
       spcnts{1,k}(i,:) = spcnts{1,k}(i,:) - gtmean(i);
   end
   %**********************************
end

%********** now show plots again, are slow trends removed???
%********** if so the units should be uncorrelated
subplot(3,3,4);
imagesc(spcnts{1,1});
title('raw cnts unit 1');
subplot(3,3,5);
imagesc(spcnts{1,2});
title('raw cnts unit 2');
subplot(3,3,6);
s1 = reshape(spcnts{1,1},[1 (T*B)]);
s2 = reshape(spcnts{1,2},[1 (T*B)]);
[r,p] = corrcoef(s1,s2);
plot(s1,s2,'ko');
title(sprintf('Raw Corr r=%+5.3f (p=%6.4f)',r(1,2),p(1,2)));



return;

