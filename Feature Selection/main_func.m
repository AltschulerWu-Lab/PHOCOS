close all
clear all

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


figure('Name', 'Figure 4C')
KL=kmeans(DataMat,5,'replicates',200);
PlotEntropyMap(DataMat',KL)

figure('Name', 'Figure 4B')
 [AIC,BIC]=CalAICBIC(DataMat);

plot([2:length(AIC)],(log10(AIC(2:end))),'b-')
hold on
plot([2:length(BIC)],(log10(BIC(2:end))),'r-')

ylabel('log10(Indicator)')
xlabel('number of classes')

legend('AIC','BIC')


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

%Plot_Label=1 : Brightest pixel class (B)
%Plot_Label=2 : Intensity class (I)
%Plot_Label=3 : Morphology class (M)
%Plot_Label=4 : Polarity class (P)
%Plot_Label=5 : Texture class (T)


for i=1:1:max(KL)
pos=find(a1==i)';
if pos(1)==1
Plot_Label(pos)=2;
elseif pos(1)==9
Plot_Label(pos)=3;
elseif pos(1)==20
Plot_Label(pos)=4;
elseif pos(1)==43
Plot_Label(pos)=5;
elseif pos(1)==15
Plot_Label(pos)=1;
end
end








for i=1:1:max(KL)
    pos=find(Plot_Label==i)';

current_attribute=Attribute(pos);
[~,sortIdx]=sort(current_attribute);

newpos=pos(sortIdx);

 

if i==1
    Image=(GR(pos,:)');
    Order=pos;
else
    Image=[Image,(GR(pos,:)')];
    Order=[Order;pos];


end

end

figure('Name', 'Figure 4D')

subplot(1,2,1)
imagesc(Image')
[a,b]=size(Image);
colormap('gray')
for i=1:1:b
cmbFeatureName{i}=[FeatureName{Order(i)}];
cmbL(i,:)=Label(Order(i),:);
end

set(gca,'ytick',1:b)
set(gca,'YTickLabel',cmbFeatureName)



set(gca,'FontSize',12);
subplot(1,2,2)
imagesc(cmbL)
caxis([-4,4])
colormap('gray')

save([intermediate_result_file,'Graph Inference\','data_five_features'],'X_LatA','X_Jas','X_Noco','X_Taxol','X_Y27','X_Calp')

Representative_Feature_Name=FeatureName(SelPos)




