clear all
close all
lambda=0.5;

for ss=1:1:7
for ii=1:1:10
for jj=1:1:10
SparseRatio=0.55+0.05*ss;

MissingRatio=0.02*ii;   
NoiseRatio=0.03*jj;

nodes_num=3;
feature_per_layer=5;

N=15;
         





for cnt=1:1:10

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

for marker=1:1:nodes_num
    
    for iii=(marker-1)*feature_per_layer+1:1:(feature_per_layer)*marker
        for jjj=1:1:feature_per_layer
            i_cnt=iii;
            j_cnt=(marker-1)*feature_per_layer+jjj;
    if G_sparse(i_cnt,j_cnt)~=0 
    G_sparse(i_cnt,j_cnt)=0;
    end
        end
    end

end




G_sparse(find(G<1e-3))=0;%
      

[s,v,d]=svd(G_sparse);
if max(max(v))>1
G_sparse=G_sparse./(max(max(v))+0.1);
end
I=eye(size(G_sparse));
D_cmb=G_sparse*pinv(I-G_sparse);% D_theo is the theoretic obeserved graph


% Random choose Missing Positions
D_vec=reshape(D_cmb,1,numel(D_cmb));
randpos=randperm(numel(D_vec));
endpos=floor(numel(D_vec)*NoiseRatio);% 
Noisepos=randpos(1:endpos); 

randN=0.2*randn(numel(D_cmb),1);
Noise=zeros(numel(D_cmb),1);
Noise(Noisepos)=randN(Noisepos);
D_vec=D_vec+Noise';
D_vec(find(D_vec<0))=0;
D=reshape(D_vec,N,N);
D(find(D<0.1*mean(mean(D))))=0;



MissingPos=[];
M_mat=zeros(size(D));
for marker=1:1:nodes_num
    
    for iii=(marker-1)*feature_per_layer+1:1:(feature_per_layer)*marker
        for jjj=1:1:feature_per_layer
            i_cnt=iii;
            j_cnt=(marker-1)*feature_per_layer+jjj;
    if D(i_cnt,j_cnt)~=0 
    M_mat(i_cnt,j_cnt)=-D(i_cnt,j_cnt);
    end
        end
    end

end



D=D+M_mat;% remove intra-marker links


D_vec=reshape(D,numel(D),1);% Reformulate the G as a vector

link_pos=find(D_vec>0);
randpos=randperm(numel(link_pos));
endpos=floor(numel(find(D_vec>0))*MissingRatio);% 
Missing_pos=link_pos(randpos(1:endpos)); 

M_vec2=zeros(size(D_vec));
M_vec2(Missing_pos)=-D_vec(Missing_pos);

M_mat2=reshape(M_vec2,N,N);

D=D+M_mat2;



D(find(D<0))=0;







D=D.*(ones(size(D))-I);




[G_r,E_r]=Robust_Direct_Effect_Pursuit(D,30);
[G_r_c]=ClosedFormEstimation(D);
[~,v,~]=svd(D);

thre_G=mean(G_r(find(G_r>0)));
thre_G_c=mean(G_r_c(find(G_r_c>0)));


G_r(find(G_r<0.5*thre_G))=0;
G_r_c(find(G_r_c<0.5*thre_G_c))=0;


[TP,TN,FN,FP]=Evaluation(G_r,G_sparse);
[TP2,TN2,FN2,FP2]=Evaluation(G_r_c,G_sparse);

[P1(cnt),R1(cnt),F1(cnt)]   =CalculatePecision(TP,FP,TN,FN);

[P2(cnt),R2(cnt),F2(cnt)]=CalculatePecision(TP2,FP2,TN2,FN2);


end
Result1(ii,jj).P=[mean(P1),std(P1)];
Result1(ii,jj).R=[mean(R1),std(R1)];
Result1(ii,jj).F=[mean(F1),std(F1)];
Result1(ii,jj).SparseRatio=SparseRatio;

Result2(ii,jj).P=[mean(P2),std(P2)];
Result2(ii,jj).R=[mean(R2),std(R2)];
Result2(ii,jj).F=[mean(F2),std(F2)];
Result2(ii,jj).SparseRatio=SparseRatio;


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


for ii=1:1:10
    for jj=1:1:10
        datamat2(10-ii+1,jj)=Result2(ii,jj).F(1);
    end
end
Acc1(ss)=mean(mean(datamat));
std1(ss)=std(reshape(datamat,numel(datamat),1));
Acc2(ss)=mean(mean(datamat2));
std2(ss)=std(reshape(datamat2,numel(datamat2),1));
end

for ii=1:1:7
label_x{ii}=num2str(0.55+0.05*ii);
end
figure('Name', 'Figure 4E')
fontsize=40;
plot(Acc1,'r-','linewidth',4)
hold on
plot(Acc2,'b-','linewidth',4)
legend('PHOCOS','CF')
 set(gca,'XTickLabel',label_x,'FontSize',fontsize);

