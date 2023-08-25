function [TAC_i,S_i]= BCoSaMP_common(H,M,Y,nTx_AS,nTx,nRx,N0,timeslot) 

TAC_i=[];
V_s=[];
R_i=Y;
K=nTx_AS;

D=kron(H,eye(timeslot));  
p=reshape(Y.',[],1);


for i=1:1:nTx_AS
    
    R_test=H'*R_i;
    R_test=R_test.';
    Dist_R=[];

    for j=1:1:nTx
        B=norm(R_test(:,j),'fro')^2;
        Dist_R=[Dist_R,B];  
    end
    
    %选出内积值最大的K列
    [V,l]=sort(Dist_R,'descend');
    diff_i=setdiff(l,TAC_i,'stable');
    TAC_Can= l(1:2*K);                 
    TAC_i = union(TAC_i,TAC_Can);   
    
    V_s=[V_s,V(1)];
    [TAC_i,v]=sort(TAC_i,'ascend');  
    
    Dist_S=[];
    if length(TAC_i)<=nTx
           H_i=H(:,TAC_i);
    end
        W=pinv(H(:,TAC_i));
        H_i=H(:,TAC_i);
        S_i=W*Y;                         
        S_re=S_i.';
        S_all=zeros(timeslot,nTx);
        S_all(:,TAC_i)=S_re;

        for  j=1:1:2*nTx_AS
            B_S=norm(S_re(:,j),'fro')^2;
            Dist_S=[Dist_S,B_S];  
        end
           [V_Si,l_Si]=sort(Dist_S,'descend');
           TAC_i=TAC_i(l_Si(1:K));
           [TAC_i,v]=sort(TAC_i,'ascend');   
           
           H_i=H(:,TAC_i);
           S_i_f=S_all(:,TAC_i);
           
           S_i_de=reshape(S_i_f,[],1);
           D_mod=rx_demodulation(S_i_de,M(1));
           S_i=tx_modulation(D_mod,M(1));
           S_i=reshape(S_i,timeslot,[]);
           S_i=S_i.';
           R_i=R_i-H_i*S_i;
end

