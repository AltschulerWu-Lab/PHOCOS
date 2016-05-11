clear all
close all

intermediate_result_file='W:\2015_01_Phenotypic_Crosstalk\code for publication\Intermediate results\'
load ([intermediate_result_file,'Graph Inference\','graph_result'])

G_i=G{1};
Gi_1=CalPersistentGraph(G_i.TG(1));
Gi_2=CalPersistentGraph(G_i.TG(2));
Gi_3=CalPersistentGraph(G_i.TG(3));
Gi_4=CalPersistentGraph(G_i.TG(4));
Gi_5=CalPersistentGraph(G_i.TG(5));

featurenum=size(Gi_1,1)/3;



 R1=CalLinkPatterns(Gi_1,featurenum);
 R2=CalLinkPatterns(Gi_2,featurenum);
 R3=CalLinkPatterns(Gi_3,featurenum);
R4=CalLinkPatterns(Gi_4,featurenum);
 R5=CalLinkPatterns(Gi_5,featurenum);
 
 
 N_R1=R1./sum(sum(R1));
 N_R2=R2./sum(sum(R2));
 N_R3=R3./sum(sum(R3));
 N_R4=R4./sum(sum(R4));
 N_R5=R5./sum(sum(R5));
 
  % The weigting matrix record the markers in the order of F,M,B
 N_R1(find(N_R1<=0.2))=0;
 N_R2(find(N_R2<=0.2))=0;
 N_R3(find(N_R3<=0.2))=0;
 N_R4(find(N_R4<=0.2))=0;
 N_R5(find(N_R5<=0.2))=0;