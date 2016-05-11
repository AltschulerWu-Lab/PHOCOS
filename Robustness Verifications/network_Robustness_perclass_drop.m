clear all
GenerateKnearestCenters()
currentFolder = pwd;
pos=find(currentFolder==filesep);
currentFolder(pos(end)+1:end)=[];
pos=[];
intermediate_result_file=[currentFolder,'Intermediate results',filesep];


for ii=1:5% select which feature class to replace
    
    for jj=1:1:5 % select which distant feature to use
        
read_name=strcat([intermediate_result_file,'Robustness Verification',filesep],'data_center_',num2str(jj));
load(read_name)

tmp_X_LatA=X_LatA;
tmp_X_Jas=X_Jas;
tmp_X_Noco=X_Noco;
tmp_X_Taxol=X_Taxol;
tmp_X_Calp=X_Calp;
tmp_X_Y27=X_Y27;

load([intermediate_result_file,'Robustness Verification',filesep,'data_center_1']);

X_LatA.ACT(ii,:)=tmp_X_LatA.ACT(ii,:);
X_LatA.MT(ii,:)=tmp_X_LatA.MT(ii,:);
X_LatA.PMY(ii,:)=tmp_X_LatA.PMY(ii,:);

X_Jas.ACT(ii,:)=tmp_X_Jas.ACT(ii,:);
X_Jas.MT(ii,:)=tmp_X_Jas.MT(ii,:);
X_Jas.PMY(ii,:)=tmp_X_Jas.PMY(ii,:);



X_Noco.ACT(ii,:)=tmp_X_Noco.ACT(ii,:);
X_Noco.MT(ii,:)=tmp_X_Noco.MT(ii,:);
X_Noco.PMY(ii,:)=tmp_X_Noco.PMY(ii,:);

X_Taxol.ACT(ii,:)=tmp_X_Taxol.ACT(ii,:);
X_Taxol.MT(ii,:)=tmp_X_Taxol.MT(ii,:);
X_Taxol.PMY(ii,:)=tmp_X_Taxol.PMY(ii,:);

X_Calp.ACT(ii,:)=tmp_X_Calp.ACT(ii,:);
X_Calp.MT(ii,:)=tmp_X_Calp.MT(ii,:);
X_Calp.PMY(ii,:)=tmp_X_Calp.PMY(ii,:);

X_Y27.ACT(ii,:)=tmp_X_Y27.ACT(ii,:);
X_Y27.MT(ii,:)=tmp_X_Y27.MT(ii,:);
X_Y27.PMY(ii,:)=tmp_X_Y27.PMY(ii,:);


beta=100;

sel_fea_class=[1:5];
n=size(sel_fea_class,2);% Node number per marker

for t=1:1:5
%%%%%%Write drug effects (LatA and Jas) into matrix    


NodeMat=[X_LatA.ACT(sel_fea_class,t);X_LatA.MT(sel_fea_class,t);X_LatA.PMY(sel_fea_class,t)];
drug_dir=zeros(size(NodeMat));
NodeMat=[NodeMat,[X_Jas.ACT(sel_fea_class,t);X_Jas.MT(sel_fea_class,t);X_Jas.PMY(sel_fea_class,t)]];
drug_dir=zeros(size(NodeMat,1),1);
drug_dir(1:n)=1;
DrugTarget=repmat(drug_dir,1,2);


%%%%%%Write drug effects (Noco and Taxol) into matrix    

NodeMat=[NodeMat,[X_Noco.ACT(sel_fea_class,t);X_Noco.MT(sel_fea_class,t);X_Noco.PMY(sel_fea_class,t)]];
drug_dir=zeros(size(NodeMat));
NodeMat=[NodeMat,[X_Taxol.ACT(sel_fea_class,t);X_Taxol.MT(sel_fea_class,t);X_Taxol.PMY(sel_fea_class,t)]];

 drug_dir=zeros(size(NodeMat,1),1);
drug_dir(n+1:2*n)=1;
DrugTarget=[DrugTarget,repmat(drug_dir,1,2)];



%%%%%%Write drug effects (Calp and Y27) into matrix    

NodeMat=[NodeMat,[X_Calp.ACT(sel_fea_class,t);X_Calp.MT(sel_fea_class,t);X_Calp.PMY(sel_fea_class,t)]];
drug_dir=zeros(size(NodeMat));
NodeMat=[NodeMat,[X_Y27.ACT(sel_fea_class,t);X_Y27.MT(sel_fea_class,t);X_Y27.PMY(sel_fea_class,t)]];
drug_dir=zeros(size(NodeMat,1),1);
drug_dir(2*n+1:3*n)=1;
DrugTarget=[DrugTarget,repmat(drug_dir,1,2)];

X{t}.data=NodeMat;
X{t}.DrugTarget=DrugTarget;
    
end

[D,D_p]=Infer_Cross_Phenotype_Graph(X);
GG{1}=D;






drawflag=false;
[G]=PHOCOS_Graph_Reduction(GG,beta,drawflag)

G_sparse=G{1}.PG;

savefile=[intermediate_result_file,'Robustness Verification',filesep,filesep,'cmpset_',num2str(ii),'_',num2str(jj),'.mat'];
save(savefile,'G','D');
    end
end

for ii=1:1:5 % which class to replace
for jj=1:1:5 % which distant to use
cur_graph=[intermediate_result_file,'Robustness Verification',filesep,'cmpset_',num2str(ii),'_',num2str(jj)];
load(cur_graph)
D1=D;
G1=G;

load([intermediate_result_file,'Graph Inference',filesep,'graph_result']); % Load the refered graph for comparisons

for i=1:1:5

    
    [A(i),B(i),C(i)]=EvaluateOverLap(D{i},D1{i});
    Indicator(i)=B(i)/(A(i)+B(i)+C(i)+B(i)-B(i));
    [A2(i),B2(i),C2(i)]=EvaluateOverLap(G{1}.TG{i},G1{1}.TG{i});
    Indicator2(i)=B2(i)/(A2(i)+B2(i)+C2(i)+B2(i)-B2(i));


end
mean_indicator_D(ii,jj)=mean(Indicator);
mean_indicator_G(ii,jj)=mean(Indicator2);
end
end


figure('Name', 'Figure 7B Raw')
plot(mean_indicator_D','linewidth',12,'marker','o','markersize',12)
legend('I','M','B','P','T')
set(gca,'fontsize',42)
axis([1,5,0.5,1])
figure('Name', 'Figure 7B PHOCOS')
plot(mean_indicator_G','linewidth',12,'marker','o','markersize',12)
legend('I','M','B','P','T')
set(gca,'fontsize',42)
axis([1,5,0.5,1])





