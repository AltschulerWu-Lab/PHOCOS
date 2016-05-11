function [a,b,c]=EvaluateOverLap(G_ref,G_sparse)

if size(G_sparse,1)==12
G_ref([5,10,15],:)=[];
G_ref(:,[5,10,15])=[];

end
G_vec1=reshape(G_ref,numel(G_ref),1);
G_vec2=reshape(G_sparse,numel(G_sparse),1);
[G_common]=intersect(find(G_vec1>0),find(G_vec2>0));
b=length(G_common);

a=sum(sum(abs(sign(G_vec1))))-b;
c=sum(sum(abs(sign(G_vec2))))-b;





