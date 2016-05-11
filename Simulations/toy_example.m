  clear all
close all

%If you want to visualize the graph drawing results, please
%install the "graphViz4Matlab" function into your local machine via:
%http://www.mathworks.com/matlabcentral/fileexchange/21652-graphviz4matlab
%Then, you can set drawflag to true for visualization purpose.
% We default set the drawflag to false that allows you running our codes
%without graph plotting function.
%The default condition will only return the link matrix and its corresponding nodes names.

drawflag=false;


feature_per_layer=5;
nodes_num=3;

N=feature_per_layer*nodes_num;



SparseRatio=0.88;
NoiseRatio=0.2;%Noise strength;
MissingLinkRatio=0.1;%Missing link ratio strength;

rng(1)
G=1*rand(N,N);
        

%Generate a sparse graph with SparseRatio
for i=1:1:size(G,1)
G(i,i)=0;
end

G=reshape(G,numel(G),1);% Reformulate the G as a vector
rng(2)
randpos=randperm(numel(G));
endpos=floor(numel(G)*SparseRatio);
removepos=randpos(1:endpos);
G(removepos)=0;
G_sparse=reshape(G,N,N);% The ground truth Graph



      

Pos_Mat=zeros(size(G_sparse));% The paiwise link positons in the same biomarker 
for marker=1:1:nodes_num
    
    for ii=(marker-1)*feature_per_layer+1:1:(feature_per_layer)*marker
        for jj=1:1:feature_per_layer
            i_cnt=ii;
            j_cnt=(marker-1)*feature_per_layer+jj;
    Pos_Mat(i_cnt,j_cnt)=1;            
        end
    end

end

G_sparse=G_sparse.*(1-Pos_Mat);

I=eye(size(G_sparse));
D_cmb=G_sparse*pinv(I-G_sparse);% D_theo is the theoretic obeserved graph


D_vec=reshape(D_cmb,numel(D_cmb),1);% Reformulate the G as a vector
rng(10)
randpos=randperm(numel(D_vec));
endpos=floor(numel(D_vec)*NoiseRatio);% 
Noisepos=randpos(1:endpos); 

rng(3)
randN=0.2*randn(numel(D_cmb),1);
Noise=zeros(numel(D_cmb),1);
Noise(Noisepos)=randN(Noisepos);
D_vec=D_vec+Noise;
D_vec(find(D_vec<0))=0;
D=reshape(D_vec,N,N);
D(find(D<0.1*mean(mean(D))))=0;

% Randomly choose Missing Positions







MissingPos=[];
M_mat=zeros(size(D));
for marker=1:1:nodes_num
    
    for ii=(marker-1)*feature_per_layer+1:1:(feature_per_layer)*marker
        for jj=1:1:feature_per_layer
            i_cnt=ii;
            j_cnt=(marker-1)*feature_per_layer+jj;
    if D(i_cnt,j_cnt)~=0 
    M_mat(i_cnt,j_cnt)=-D(i_cnt,j_cnt);
    end
        end
    end

end



removed_missing_ratio=numel(find(M_mat)~=0)/numel(find(D_cmb)~=0);
remaining_missin_ratio=max(0,MissingLinkRatio-removed_missing_ratio);

D=D+M_mat;

D_vec=reshape(D,numel(D),1);% Reformulate the G as a vector
rng(4)
link_pos=find(D_vec>0);
randpos=randperm(numel(link_pos));
endpos=floor(numel(find(D_vec>0))*remaining_missin_ratio);% 
Missing_pos=link_pos(randpos(1:endpos)); 

M_vec2=zeros(size(D_vec));
M_vec2(Missing_pos)=-D_vec(Missing_pos);

M_mat2=reshape(M_vec2,N,N);

D=D+M_mat2;

E=M_mat+M_mat2;

pos=find(D>0);


D=D.*(ones(size(D))-I);


[G_r,E_r]=Robust_Direct_Effect_Pursuit(D,50);

pos=find(D==0);
G_r(pos)=0;
[G_r_c]=ClosedFormEstimation(D);
[~,v,~]=svd(D);

pos=find(G_r>0);
meanval=mean(mean(G_r(pos)));

G_r(find(G_r)<meanval)=0;

G_sparse=G_sparse./max(max(G_sparse));

G_r=G_r./max(max(G_r));


thre_G=mean(G_r(find(G_r>0)));
thre_G_c=mean(G_r_c(find(G_r_c>0)));


G_r(find(G_r<1*thre_G))=0;
G_r_c(find(G_r_c<0.8*thre_G_c))=0;
E_r(find(G_r>0))=0;
E_r=abs(E_r);
thre_E=mean(E_r(find(E_r>0)));

 E_r(find(abs(E_r)<0.8*thre_E))=0;
 
 
 
E=abs(E)./max(max(abs(E)));

E_r=E_r./max(max(E_r));


MissingLinks=E_r.*Pos_Mat;
val=sort(reshape(MissingLinks,1,numel(MissingLinks)),'descend');



DrawGraph_toy(G_sparse,drawflag) % Ground truth
 
DrawGraph_toy(D,drawflag) % Observed graph

DrawGraph_toy(G_r_c,drawflag) % CF recovery

DrawGraph_toy(G_r,drawflag) % PHOCOS recovery


[TP,TN,FN,FP]=Evaluation(G_r,G_sparse);
[P1,R1,F1]=CalculatePecision(TP,FP,TN,FN);

[TP,TN,FN,FP]=Evaluation(G_r_c,G_sparse);
[P2,R2,F2]=CalculatePecision(TP,FP,TN,FN);





