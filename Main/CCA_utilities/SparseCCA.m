%%% Calculating regularized CCA, maximizing (X.wxMat)' . (Y.wyMat)
%%% Inputs X= n x p1   Y=n x p2
%%% reg -> type of regularization 1 : L1  and 2 : L2 
%%% L1 regularization parameters  cx > |wx|  and cy > |wy|  (for L2 cx and cy are both 1)
%%% numCC -> the Number of CCA Modes
%%% Output wxMat=p1 x numCC  , wyMat=p2 x numCC, rVec=numCC x 1 -> the correlation coefficient at each mode
%%% This code was developed based the the algorithm introduced by D. Witten and R. Tibshirani in the following paper:
%%% Witten, D. M. & Tibshirani, R. J. Extensions of sparse canonical correlation analysis with
%%% applications to genomic data. Stat Appl Genet Mol Biol 8, Article28, doi:10.2202/1544-6115.1470 (2009).
 
function [wxMat,wyMat,rVec]=SparseCCA(X,Y,cx,cy,reg,numCC)
dx=size(X);
dy=size(Y);
wxMat=zeros(dx(2),numCC);
wyMat=zeros(dy(2),numCC);
rVec=zeros(numCC,1);



xmask=max(X)>min(X);
ymask=max(Y)>min(Y);

X=X(:,xmask);
Y=Y(:,ymask);

X=X-repmat(mean(X),[dx(1),1]);
X=X./repmat(sqrt(var(X)),[dx(1),1]);
Y=Y-repmat(mean(Y),[dy(1),1]);
Y=Y./repmat(sqrt(var(Y)),[dy(1),1]);
tol=1e-3;
X(isnan(X))=0;
Y(isnan(Y))=0;


dx=size(X);
dy=size(Y);






for icc=1:numCC

wx=randn(dx(2),1);
wy=randn(dy(2),1);
wx=wx/norm(wx);
wy=wy/norm(wy);

wxp=zeros(dx(2),1);
wyp=zeros(dy(2),1);
iter=0;

while (norm(wx-wxp)/norm(wx)) > tol && (norm(wy-wyp)/norm(wy)) > tol
    wxp=wx;
    wyp=wy;

    wx = X' * Y *wy;
    if reg==1
        wx=softTreshold(wx,cx);
    elseif reg==2
        wx=wx/norm(wx);
    end
    
    wy=(wx' * X' * Y)';
    wy(isnan(wy))=0;
    if reg==1
        wy=softTreshold(wy,cy);
    elseif reg==2
        wy=wy/norm(wy);
    end

    iter=iter+1;
    if iter>500
        break;
    end

end
c=corrcoef([X * wx , Y * wy]);
r=c(1,2);

wxMat(xmask,icc)=wx;
wyMat(ymask,icc)=wy;
rVec(icc)=r;

X=X-(X*wx*wx');
Y=Y-(Y*wy*wy');

end  


end