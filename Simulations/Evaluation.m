function [TP,TN,FN,FP]=Evaluation(G_infer,G_gnd)

DiffMat=sign(G_infer)-sign(G_gnd);
FP=length(find(DiffMat==1));
FN=length(find(DiffMat==-1));
pos1=find(G_gnd==0);
pos2=find(G_infer==0);

TN=length(intersect(pos1,pos2));

pos3=find(sign(G_gnd)==1);
pos4=find(sign(G_infer)==1);

TP=length(intersect(pos3,pos4));

