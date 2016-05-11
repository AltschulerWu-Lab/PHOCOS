function G=ClosedFormEstimation(D)
max_eig=1;
[~,v,~]=svd(D);
if v(1,1)>max_eig
coef=v(1,1)/max_eig;
D=D/coef;
end

I=eye(size(D));
G=pinv(I+D)*D;