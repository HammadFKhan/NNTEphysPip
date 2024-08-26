%%%%Simulations Outputs
TransMtx={};
AdjMtx={};
CCAModesX={};CCAModesY={};
PCAModes={};
CCAVal={};
CCATrain={};
% %%%%%%%%%%%Iterations
% CoupVec=[0.26,0.15,0.11,0.085,0.08];
Iter=0;

for K=[2,4,6,8,10]
    Iter=Iter+1
    Iterb=0;
    for beta=[0,0.2,0.4,0.6,0.8,1]
        Iterb=Iterb+1
        x={};
        xn={};
        W={};
        N=11;
        cellPerArea=500;
        x_int=0.001/5;
        w_int=0.9/sqrt(K);%CoupVec(Iter);
        sOut=0.01/5;%%% internal random added noise
        TMax=50000;
        xRecording=zeros(N,cellPerArea,TMax);
        %%%% Small-World Network Construction
        A=zeros(N);
        for i=1:N
            for j=1:N
                if (min(abs(i-j),N-abs(i-j))<=(K/2) && i~=j)
                    A(i,j)=1;
                end
            end
        end
        for i=1:N
            for j=1:N
                if (A(i,j)==1 && rand(1)<beta  && i~=j)
                    A(i,j)=0;
                    A(i,i)=1;
                    OpenNodes=find(A(i,:)==0);
                    A(i,i)=0;
                    newNode=randperm(length(OpenNodes),1);
                    A(i,OpenNodes(newNode))=1;
                    
                end
            end
        end
        A=A';
        AdjMtx{Iter,Iterb}=A;
        %%%%Initialization
        for i=1:N
            x{i}=x_int * randn(cellPerArea,1);
            for j=1:N
                if A(i,j)==1
                    W{i,j}=randn(cellPerArea,cellPerArea);
                    for nc=1:cellPerArea
                        W{i,j}(nc,:)=w_int*( W{i,j}(nc,:) / norm(W{i,j}(nc,:)));
                    end
                end
                
            end
        end
         TransMtx{Iter,Iterb}=W;
%         %%%%Network Stimulation
        for tc=1:TMax
            for j=1:N
                xRecording(j,:,tc)=x{j};
                
                xn{j}=(sOut*randn(cellPerArea,1)/10) ;
                for i=1:N
                    if A(i,j)==1
                        xn{j}=xn{j}+W{i,j}*x{i};
                    end
                end
                %xn{j}=xn{j}- (w_int*sum(A(:,j))*x{j});

            end
            x=xn;
        end
        assert(max(max(max(xRecording)))<1000,'Network is not stable');
%         %%%Fluctuation Analysis
        Wx={};
        Wy={};
        rTrain=zeros(N,N,20);
        rVal=zeros(N,N,20);
        
        for i=1:N
            i;
            for j=1:N
                j;
                if i<j
                    Wx{i,j}=zeros(cellPerArea,20);
                    Wy{i,j}=zeros(cellPerArea,20);
                    [Wx{i,j},Wy{i,j},rTrain(i,j,:)]=SparseCCA(squeeze(xRecording(i,:,500:25000))',squeeze(xRecording(j,:,500:25000))',1,1,2,20);
                    for cc=1:20
                        t1=squeeze(xRecording(i,:,25000:50000))'*squeeze(Wx{i,j}(:,cc));
                        t2=squeeze(xRecording(j,:,25000:50000))'*squeeze(Wy{i,j}(:,cc));
                        c = abs(corrcoef([t1,t2]));
                        r=c(1,2);
                        rVal(i,j,cc)=r;
                    end
                    
                end
            end
        end
        CCAModesX{Iter,Iterb}=Wx;
        CCAModesY{Iter,Iterb}=Wy;
        
        CCAVal{Iter,Iterb}=rVal;
        CCATrain{Iter,Iterb}=rTrain;
        %%%%%%%%%%%%%%%%% PCA Analysis
        
        [COEFF, SCORE, LATENT] = pca(reshape(xRecording(:,:,500:end),[N*cellPerArea,TMax-499])','NumComponents',20);
        
        PCAModes{Iter,Iterb}=COEFF;
    end
end
%%%%%%%%%%Save Results
 save(FileName,'sOut','x_int','TransMtx','AdjMtx','CCAModesX','CCAModesY','PCAModes','CCAVal','CCATrain','-v7.3');

%%%%%%%%%%%%%%%% Visualization
%% CCAVal
cc=1;
res=[];
for k=1:5
    buf=[];
    for b=1:6
        buf=[buf;zeros(1,N);(squeeze(CCAVal{k,b}(:,:,cc))+squeeze(CCAVal{k,b}(:,:,cc))')];
    end
    res=[res,zeros(72,1),buf];
end
figure();colormap('jet');imagesc(res);
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
    set(gca,'ytick',[])
set(gca,'yticklabel',[])
%%%%%%%%%%%%%%
%%%%CCA Modes X,Y
cc=1;
SimMtx=zeros(N,20,6,N,N,5);
for I1=1:5
    I1
    for I2=1:6
        I2
Wx=CCAModesX{I1,I2};
Wy=CCAModesY{I1,I2};

for ccn=1:20
for s=1:N
for d1=1:N
    for d2=1:N
        if s~=d1 && s~=d2 && d1~=d2
            if s<d1
                w1=Wx{s,d1}(:,ccn);
            else
                w1=Wy{d1,s}(:,ccn);
            end
            if s<d2
                w2=Wx{s,d2}(:,ccn);
            else
                w2=Wy{d2,s}(:,ccn);
            end
            %SimMtx(d1,ccn,I2,d2,s,I1)=abs(similarity([w1,w2]));
            c = abs(corrcoef([w1,w2]));
            r=c(1,2);
            SimMtx(d1,ccn,I2,d2,s,I1)=r;
        end
    end
end
end
end
    end
end
%SimMtx=SimMtx([3,1,2,4:11],cc,:,[3,1,2,4:11],:,:);
figure();colormap('jet')
imagesc(squeeze(reshape(mean(SimMtx(:,cc,:,:,1,:),5),[N*6,N*5])))
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
    set(gca,'ytick',[])
set(gca,'yticklabel',[])

%%%%%%%%%%%%%% PCA Modes
figure();
indarea=reshape(repmat([1:20],[500,1]),[],1);
n=1;
for I1=1:5
    for I2=1:6
        PC=PCAModes{I1,I2}(:,1);
        PCMask=abs(PC-mean(PC)) > (2*sqrt(var(PC)));
        subplot(5,6,n);
        hist(indarea(PCMask),N);
        axis([1 N 0 50])
       
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
    set(gca,'ytick',[])
set(gca,'yticklabel',[])
        
        n=n+1;
    end
end
