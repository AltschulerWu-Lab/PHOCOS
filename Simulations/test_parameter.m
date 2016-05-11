clear all
close all

for tt=1:1:10
for ii=1:1:10
for jj=1:1:10
SparseRatio=0.8;

MissingRatio=0.01*ii;   
NoiseRatio=0.05*jj;
N=15;
         





for cnt=1:1:5

    G=1*rand(N,N);



G=1*rand(N,N);


%Generate a sparse graph with SparseRatio
for i=1:1:size(G,1)
G(i,i)=0;
end

G=reshape(G,numel(G),1);% Reformulate the G as a vector
randpos=randperm(numel(G));
endpos=floor(numel(G)*SparseRatio);
removepos=randpos(1:endpos);
G(removepos)=0;
G_sparse=reshape(G,N,N);% The ground truth Graph



G_sparse(find(G<1e-3))=0;%
      

[s,v,d]=svd(G_sparse);
if max(max(v))>1
G_sparse=G_sparse./(max(max(v))+0.1);
end      
I=eye(size(G_sparse));
D_cmb=G_sparse*pinv(I-G_sparse);% D_theo is the theoretic obeserved graph


% Random choose Missing Positions
D_vec=reshape(D_cmb,1,numel(D_cmb));
M_vec=zeros(size(D_vec));
N_vec=zeros(size(D_vec));



N_vec=NoiseRatio*mean(mean(D_vec))*randn(1,numel(D_vec));
N_mat=reshape(N_vec,N,N);

link_position=find(D_vec>0);
randpos=randperm(numel(link_position));
MissingNum=floor(MissingRatio*numel(randpos));
MissingPos=randpos(1:MissingNum);

M_vec(MissingPos)=-D_vec(MissingPos);
M_mat=reshape(M_vec,N,N);


E=N_mat+M_mat;


D=D_cmb+E;

D(find(D<0))=0;







D=D.*(ones(size(D))-I);



beta=10+10*tt;
[G_r,E_r]=Robust_Direct_Effect_Pursuit(D,beta);
[~,v,~]=svd(D);

thre_G=mean(G_r(find(G_r>0)));


G_r(find(G_r<0.6*thre_G))=0;


[TP,TN,FN,FP]=Evaluation(G_r,G_sparse);

[P1(cnt),R1(cnt),F1(cnt)]   =CalculatePecision(TP,FP,TN,FN);



end
Result1(ii,jj).P=[mean(P1),std(P1)];
Result1(ii,jj).R=[mean(R1),std(R1)];
Result1(ii,jj).F=[mean(F1),std(F1)];
Result1(ii,jj).SparseRatio=SparseRatio;



end


end


for ii=1:1:10
    for jj=1:1:10
        datamat(10-ii+1,jj)=Result1(ii,jj).F(1);
    end
end
cmin=2;
cmax=4;
caxis([cmin cmax])




cmin=0.5;
cmax=max([max(max(datamat))]);
fontsize=25;
caxis([cmin cmax])
for ii=1:1:10
label_x{ii}=num2str(0.05*ii);
label_y{ii}=num2str(0.55-0.05*ii);

end

mean_R(tt)=mean(mean(datamat));

end
figure('Name', 'Figure 3DC')
plot(mean_R,'b-')
hold on
legend('PHOCOS')
axis([1 10 0.5 0.9]);


