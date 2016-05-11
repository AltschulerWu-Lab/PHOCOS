function Node=DrawGraph5(G_p,drawflag)

biomarker{1}='F';
biomarker{2}='M';
biomarker{3}='B';


Fea{1}='I';
Fea{2}='M';
Fea{3}='B';
Fea{4}='P';
%Fea{5}='T';

cnt=0;
N=size(G_p,1)/3;

if N==5
Fea{5}='T';
end

NodeColor=repmat([1,0.3,0.3],N,1);
NodeColor=[NodeColor;repmat([0.2,0.5,1],N,1)];
NodeColor=[NodeColor;repmat([0.3,1,0.2],N,1)];

for i=1:1:length(biomarker)
    for j=1:1:length(Fea)
        cnt=cnt+1;
        Node{cnt}=[biomarker{i},'_',Fea{j}];
        %Node{cnt}=[Fea{j}];

    end
end

if drawflag==true
nodecolor=[]
% 
%  G_p(find(G_p<0.1))=0;
G_p=G_p./(mean(max(G_p))+0.01);

for i=1:1:N
    for j=1:1:N
%         EdgeColor{(i-1)*N+j,1}=Node{i};
%         EdgeColor{(i-1)*N+j,2}=Node{j};
        EdgeColor{(i-1)*N+j,1}=i;
        EdgeColor{(i-1)*N+j,2}=j;
%         EdgeColor{(i-1)*N+j,3}=max(((1-G_p(i,j))),0)*[1,1,1];
        EdgeColor{(i-1)*N+j,3}=[1,1,1];

    end
end


graphViz4Matlab('-adjMat',G_p, '-nodeLabels',Node,'-nodeColors', NodeColor,'-edgeColors',EdgeColor)
end