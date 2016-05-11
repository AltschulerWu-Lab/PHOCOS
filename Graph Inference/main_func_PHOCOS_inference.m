close all
clear all
%If you want to visisualize the graph drawing results, please
%install the "graphViz4Matlab" function into your local machine via:
%http://www.mathworks.com/matlabcentral/fileexchange/21652-graphviz4matlab
%Then, you can set drawflag to true for visulization purpose.
%Defaulty, we set the drawflag to false that allows you running our codes
%without graph plotting function.
%The default condition will only return the link matrix and its corresponding nodes names.
drawflag=false;

currentFolder = pwd;
pos=find(currentFolder==filesep);
currentFolder(pos(end)+1:end)=[];
pos=[];
intermediate_result_file=[currentFolder,'Intermediate results',filesep];
load ([intermediate_result_file,'Graph Inference\','data_five_features'])

beta=100; % The relaxation parameter used in PHOCOS reduction

sel_fea_class=[1:5]; % Use all the five features for graph inference
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

[D,D_p]=Infer_Cross_Phenotype_Graph(X); % D records temporary influence graphs. D_p is the persistent influence graph (before PHOCOS reduction)


Node_name=DrawGraph5(D_p,drawflag);
GG{1}=D;






[G]=PHOCOS_Graph_Reduction(GG,beta,drawflag) % This function applies PHOCOS for graph complexity reduction

G_sparse=G{1}.PG;

save([intermediate_result_file,'Graph Inference\','graph_result'],'Node_name','D','G')
