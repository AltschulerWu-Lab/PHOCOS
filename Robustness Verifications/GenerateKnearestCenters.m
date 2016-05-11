function GenerateKnearestCenters()

currentFolder = pwd;
pos=find(currentFolder==filesep);
currentFolder(pos(end)+1:end)=[];
pos=[];
intermediate_result_file=[currentFolder,'Intermediate results',filesep];
load ([intermediate_result_file,'Raw Feature Data\','FeatureName'])

load ([intermediate_result_file,'Raw Feature Data\','Jas_Fea'])
for j=1:1:length(Data.ACT);    
Feature_vector(j,:)=[Data.ACT{j}(1:5),Data.MT{j}(1:5),Data.PMY{j}(1:5)];
end
DataMat=Feature_vector;
Feature_Jas=Feature_vector;
load ([intermediate_result_file,'Raw Feature Data\','LatA_Fea'])
for j=1:1:length(Data.ACT);    
Feature_vector(j,:)=[Data.ACT{j}(1:5),Data.MT{j}(1:5),Data.PMY{j}(1:5)];
end
DataMat=[DataMat,Feature_vector];
Feature_LatA=Feature_vector;

load([intermediate_result_file,'Raw Feature Data\','Noco_Fea'])

for j=1:1:length(Data.ACT);    
Feature_vector(j,:)=[Data.ACT{j}(1:5),Data.MT{j}(1:5),Data.PMY{j}(1:5)];
end
DataMat=[DataMat,Feature_vector];
Feature_Noco=Feature_vector;

load([intermediate_result_file,'Raw Feature Data\','Taxol_Fea'])

for j=1:1:length(Data.ACT);    
Feature_vector(j,:)=[Data.ACT{j}(1:5),Data.MT{j}(1:5),Data.PMY{j}(1:5)];
end
DataMat=[DataMat,Feature_vector];
Feature_Taxol=Feature_vector;

load([intermediate_result_file,'Raw Feature Data\','Y27_Fea'])

for j=1:1:length(Data.ACT);    
Feature_vector(j,:)=[Data.ACT{j}(1:5),Data.MT{j}(1:5),Data.PMY{j}(1:5)];
end
DataMat=[DataMat,Feature_vector];
Feature_Y27=Feature_vector;



load([intermediate_result_file,'Raw Feature Data\','Calp_Fea'])

for j=1:1:length(Data.ACT);    
Feature_vector(j,:)=[Data.ACT{j}(1:5),Data.MT{j}(1:5),Data.PMY{j}(1:5)];
end
DataMat=[DataMat,Feature_vector];
Feature_Calp=Feature_vector;






load([intermediate_result_file,'Raw Feature Data\','Attribute'])  % Our annotation about feature meanings


Label=zeros(size(Attribute,1),max(Attribute));
for i=1:1:length(Attribute)
Label(i,Attribute(i))=1;
end

GR=DataMat;


KL=kmeans(DataMat,5,'replicates',200);




for i=1:1:max(KL)
pos=find(KL==i);
if pos(1)==1
SelPos(1)=findNearResposeCurve(mean(DataMat(pos,:)),DataMat,1);
elseif pos(1)==9
SelPos(2)=findNearResposeCurve(mean(DataMat(pos,:)),DataMat,1);
elseif pos(1)==15
SelPos(3)=findNearResposeCurve(mean(DataMat(pos,:)),DataMat,1);
elseif pos(1)==20
SelPos(4)=findNearResposeCurve(mean(DataMat(pos,:)),DataMat,1);
elseif pos(1)==43
SelPos(5)=findNearResposeCurve(mean(DataMat(pos,:)),DataMat,1);
end
end



for i=1:1:max(KL)
    pos=SelPos(i);
X_LatA.ACT(i,:)=Feature_LatA(pos,1:5);
X_LatA.MT(i,:)=Feature_LatA(pos,6:10);
X_LatA.PMY(i,:)=Feature_LatA(pos,11:15);

X_Jas.ACT(i,:)=Feature_Jas(pos,1:5);
X_Jas.MT(i,:)=Feature_Jas(pos,6:10);
X_Jas.PMY(i,:)=Feature_Jas(pos,11:15);


X_Noco.ACT(i,:)=Feature_Noco(pos,1:5);
X_Noco.MT(i,:)=Feature_Noco(pos,6:10);
X_Noco.PMY(i,:)=Feature_Noco(pos,11:15);


X_Taxol.ACT(i,:)=Feature_Taxol(pos,1:5);
X_Taxol.MT(i,:)=Feature_Taxol(pos,6:10);
X_Taxol.PMY(i,:)=Feature_Taxol(pos,11:15);


X_Y27.ACT(i,:)=Feature_Y27(pos,1:5);
X_Y27.MT(i,:)=Feature_Y27(pos,6:10);
X_Y27.PMY(i,:)=Feature_Y27(pos,11:15);


X_Calp.ACT(i,:)=Feature_Calp(pos,1:5);
X_Calp.MT(i,:)=Feature_Calp(pos,6:10);
X_Calp.PMY(i,:)=Feature_Calp(pos,11:15);


end





a1=KL';

Plot_Label=zeros(size(a1));




% In different runs, K-means may randomly assign a class label. We
% rearrange the class label here to make it the same as our published
% results. 

%Reorder_Label=1 : Intensity class (I)
%Reorder_Label=2 : Morphology class (M)
%Reorder_Label=3 : Brightest pixel class (B)
%Reorder_Label=4 : Polarity class (P)
%Reorder_Label=5 : Texture class (T)


  for i=1:1:max(KL)
pos=find(a1==i)';
if pos(1)==1
Reorder_Label(pos)=1;
elseif pos(1)==9
Reorder_Label(pos)=2;
elseif pos(1)==20
Reorder_Label(pos)=4;
elseif pos(1)==43
Reorder_Label(pos)=5;
elseif pos(1)==15
Reorder_Label(pos)=3;
end
end

 sel_num=5;
  
for i=1:1:max(Reorder_Label)
pos=find(Reorder_Label==i);
SelPos(i,:)=findNearResposeCurve(mean(DataMat(pos,:)),DataMat,sel_num);
end
for j=1:1:sel_num
for i=1:1:max(KL)
    pos=SelPos(i,j);
X_LatA.ACT(i,:)=mean(Feature_LatA(pos,1:5),1);
X_LatA.MT(i,:)=mean(Feature_LatA(pos,6:10),1);
X_LatA.PMY(i,:)=mean(Feature_LatA(pos,11:15),1);

X_Jas.ACT(i,:)=mean(Feature_Jas(pos,1:5),1);
X_Jas.MT(i,:)=mean(Feature_Jas(pos,6:10),1);
X_Jas.PMY(i,:)=mean(Feature_Jas(pos,11:15),1);


X_Noco.ACT(i,:)=mean(Feature_Noco(pos,1:5),1);
X_Noco.MT(i,:)=mean(Feature_Noco(pos,6:10),1);
X_Noco.PMY(i,:)=mean(Feature_Noco(pos,11:15),1);


X_Taxol.ACT(i,:)=mean(Feature_Taxol(pos,1:5),1);
X_Taxol.MT(i,:)=mean(Feature_Taxol(pos,6:10),1);
X_Taxol.PMY(i,:)=mean(Feature_Taxol(pos,11:15),1);


X_Y27.ACT(i,:)=mean(Feature_Y27(pos,1:5),1);
X_Y27.MT(i,:)=mean(Feature_Y27(pos,6:10),1);
X_Y27.PMY(i,:)=mean(Feature_Y27(pos,11:15),1);


X_Calp.ACT(i,:)=mean(Feature_Calp(pos,1:5),1);
X_Calp.MT(i,:)=mean(Feature_Calp(pos,6:10),1);
X_Calp.PMY(i,:)=mean(Feature_Calp(pos,11:15),1);


end
save([intermediate_result_file,filesep,'Robustness Verification',filesep,'data_center_',num2str(j),'.mat'],'X_LatA','X_Jas','X_Noco','X_Taxol','X_Y27','X_Calp')

end

