function G_persistent=CalPersistentGraph(G)
 for ii=1:1:size(G{1})
        for jj=1:1:size(G{1})
           for t=1:1:length(G)
           vec(t)=[G{t}(ii,jj)];
           end
           if min(min(vec))==0
             G_persistent(ii,jj)=0;
%                G_persistent(ii,jj)=mean(vec);
           else
               G_persistent(ii,jj)=max(vec);
           end
     end
 end
% G_persistent(find(G_persistent<0.05))=0;
end

