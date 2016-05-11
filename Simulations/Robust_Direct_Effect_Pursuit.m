function [G_k,E_k]=Robust_Direct_Effect_Pursuit(D,beta)
max_eig=1;
[~,v,~]=svd(D);
converge_rate=1e-5;
loop=1;
max_itr=1000;
if v(1,1)>max_eig
coef=v(1,1)/max_eig;
D=D/coef;
end

I=eye(size(D));
G_k=pinv(I+D)*D;

G_hat_k=G_k;

E_k=zeros(size(D));
E_hat_k=zeros(size(D));
%  E_hat_k=D;
G_hat_k=zeros(size(D));
E_k=E_hat_k;
G_k=G_hat_k;

 

mu=1;    
Lambda1=ones(size(G_k));
Lambda2=ones(size(G_k));


k=0;
while loop
    k=k+1;
    G_hat_kp1=shrinkage(G_k-1/mu*Lambda1,mu/2);
    E_hat_kp1=shrinkage(E_k-1/mu*Lambda2,mu/2);
% 
   E_gradient=2*beta*(E_k*(I-G_k)+D*(I-G_k)-G_k)*(I-G_k)'+mu*(E_k-(E_hat_kp1+mu^-1*Lambda2));
   G_gradient=2*beta*(I+E_k+D)'*((I+E_k+D)*G_k-(D+E_k))+mu*(G_k-(G_hat_kp1+mu^-1*Lambda1));
    lr=0.001*(0.99)^k;   
   E_kp1=E_k-lr*E_gradient;
   G_kp1=G_k-lr*G_gradient;

  G_k=G_kp1;
  E_k=E_kp1;

  Lambda1=Lambda1+mu*(G_hat_kp1-G_kp1);
  Lambda2=Lambda2+mu*(E_hat_kp1-E_kp1);

  
  mu=min(1.1*mu,100000);

  E_k(find((E_k)<0))=0;
   G_k(find(G_k<0))=0;
  G_hat_kp1=G_hat_k;
  E_hat_kp1=E_hat_k;
 
 
 
  
  obj(k)=sum(sum(abs(G_k)))+sum(sum(abs(E_k)))+beta*norm((D+E_k)*(I-G_k)-G_k);
  if k>1 &&abs(obj(k)-obj(k-1))/obj(k-1)<converge_rate||k>max_itr
%       plot(obj(1:50))
      loop=0;
  end
end


  G_k=G_k.*(ones(size(G_k))-I);
  E_k=E_k.*(ones(size(E_k))-I);
    
end
    
    function B=shrinkage(A,alpha)

        pos=find(abs(A)<1./alpha);
        B=sign(A).*(abs(A)-1./alpha);
        B(pos)=0;
        

    end
    
        function obj=Calobj(B,A,E,beta)
        I=eye(size(A));
        obj=1/2*norm(B-pinv(I-A)-E)+beta*sum(sum(abs(E)));
        

    end