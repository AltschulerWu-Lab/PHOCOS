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
load ([intermediate_result_file,'Graph Inference\','graph_result']);

G_i=G{1};
 Gi_ini=CalPersistentGraph(G_i.TG(1:2));
 Gi_est=CalPersistentGraph(G_i.TG(3:4));
 Gi_mainten=CalPersistentGraph(G_i.TG(5));

 Gi_p=CalPersistentGraph(G_i.TG(1:5));

featurenum=size(Gi_ini,1)/3;

% 
 Node_Name=DrawGraph5(Gi_ini,drawflag);
 
 DrawGraph5(Gi_est,drawflag);
 
 DrawGraph5(Gi_mainten,drawflag);



 R1=CalLinkPatterns(Gi_ini,featurenum);
 R2=CalLinkPatterns(Gi_est,featurenum);
 R3=CalLinkPatterns(Gi_mainten,featurenum);
 
 N_R1=R1./sum(sum(R1));
 N_R2=R2./sum(sum(R2));
 N_R3=R3./sum(sum(R3));
 
 % The weigting matrix record the module level interactions
 % Three modules are arranged in the order of F,M,B
 % See Fig6.b in Deng et al. "PHOCOS: Inferring Multi-Feature Phenotypic
 % Crosstalk Networks" Bioinformatics, 2016
 N_R1(find(N_R1<=0.2))=0;
 N_R2(find(N_R2<=0.2))=0;
 N_R3(find(N_R3<=0.2))=0;

