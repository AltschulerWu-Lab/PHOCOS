function [D,D_persistent]=Infer_Cross_Phenotype_Graph(X)

threshold=0.25;
for t=1:1:5
Data_cur=X{t}.data;

DrugDir=X{t}.DrugTarget;

for i=1:1:size(DrugDir,2)
    cur_dir=DrugDir(:,i);
 cur_drug=abs(Data_cur(:,i));
 node_response=abs(Data_cur(:,i));
 drug_Graph{i}=zeros(size(node_response,1),size(node_response,1));
  
     max1=4;
     node_response=node_response./max1;% normalize the z-score response into the range of [0,1]
     cur_drug=cur_drug./max1;% normalize the z-score response into the range of [0,1]
     
     pos=find(cur_drug<threshold);
     cur_drug(pos)=0;% Remove the responses on the nodes whose zscores are tiny (not actived)
     pos=find(node_response<threshold);
     node_response(pos)=0;% Remove the responses on the nodes whose zscores are tiny (not actived)
     drug_node=find(cur_dir==1);
     
 for ii=1:1:length(node_response)
    
     if cur_dir(ii)~=1% not directly drug targeted 
         if abs(node_response(ii))>threshold
             
             active_pos=(drug_node);% Link should be definetely predicted from drug perturbed nodes to the current nodes
             weight=zeros(size(node_response));
             weight(active_pos)=node_response(active_pos);

             weight=node_response(ii).*weight/sum(cur_drug(active_pos));            
             drug_Graph{i}(:,ii)=weight;
             
         end
     end
     
     
 end
 
 
end

%% Merge the links found with different drug perturbations
    for ii=1:1:size(drug_Graph{1})
        for jj=1:1:size(drug_Graph{1})
           for n=1:1:length(drug_Graph)
           vec(n)=[drug_Graph{n}(ii,jj)];
           end
            Graph(ii,jj)=max(vec);
     end
    end
    


% Graph(find(Graph<0.05))=0;
D{t}=Graph;
end

    vec=zeros(1,length(D));
    for ii=1:1:size(D{1})
        for jj=1:1:size(D{1})
           for t=1:1:length(D)
           vec(t)=[D{t}(ii,jj)];
           end
  if min(min(vec))==0
             D_persistent(ii,jj)=0;
%                G_persistent(ii,jj)=mean(vec);
           else
               D_persistent(ii,jj)=mean(vec);
  end
  
        end
end
