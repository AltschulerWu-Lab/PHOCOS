function [G]=PHOCOS_Graph_Reduction(Data,beta,drawflag)

for i=1:1:length(Data)
D=Data{i};
FG_i=[];
for t=1:1:length(D)
   [T1,E]=Robust_Direct_Effect_Pursuit(D{t},beta);
 FG_i{t}=CalGraph(T1);

end




Gp_1=CalPersistentGraph(FG_i);% remove  very small weights after sparse optimization

   G{i}.TG=FG_i;
   G{i}.PG=Gp_1;
  DrawGraph5(Gp_1,drawflag);
end


% Gp_1=Gp_1./max(max(Gp_1));
% Gp_2=Gp_2./max(max(Gp_2));

% for i=1:1:size(Gp_1,1)
%     for j=1:1:size(Gp_2,1)
%     if Gp_1(i,j)>Gp_2(i,j)
%         Gp_2(i,j)=0;
%     else
%          Gp_1(i,j)=0;
%     end
%     end
% end
% Gp_1=CalStrongestInteraction(Gp_1,n);
% Gp_1=CalStrongestInteraction(Gp_1,n);

% 
%  DrawGraph5(Gp_1)
% 
%   DrawGraph5(Gp_2)


   



end



function G=CalGraph(G)
pos=find(G<0.6*mean(mean(G(find(G>0)))));
G(pos)=0;
N=size(G)/3;
for i=1:1:3
    for jj=(i-1)*N+1:1:i*N
        for kk=(i-1)*N+1:1:i*N
            G(jj,kk)=0;
        end
    end
end
end



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
end



