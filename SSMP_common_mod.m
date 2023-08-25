function  [detect_bit_SSMP,N_error,S_loc,TAC_final,S_Set]=SSMP_common_mod(H,Y,nRx,nTx,N0,n_tx,M,nTx_AS,TAC_right)

R_y=Y*Y'/7;  %为什么要除以7
[V,D]=eig(R_y);%特征值构成对角矩阵D，对应特征向量构成对角矩阵V
N_error=0;
S_Set=[];
U_y=[];

for r=nRx:-1:1
    if D(r,r)>0.01      
        U_y=[U_y,V(:,r)];  %选择支撑集
    end
end

P_i=eye(nRx);

T_set=[];
T_index=[];
R_0=Y;
A=[1:1:nTx];              

dist=zeros(1,nTx);
for i_iter=1:nTx_AS                 
    T_set_left=setdiff(A,T_set,'stable');
    
    S_left=nTx-length(T_set); %剩余支撑集的大小；
    dist(T_set)=0;
    
    for i1=1:S_left
        k=T_set_left(i1);            %T_set_left是T_set中不涉及的索引
        n_k=(k-1)*n_tx+1:k*n_tx;
        H_k_i=P_i*H(:,n_k);
        dist(k)=norm(U_y'*H_k_i)^2/(norm(H_k_i)^2);
    end
    [v,k_hat]=max(dist);             %找到剩余集中最有可能的一个
    [v_test,k_test]=sort(dist);
    T_set=[T_set,k_hat];
    
    B=H(:,T_index)*pinv(H(:,T_index));
    P_i=eye(nRx)-B;
    R_i=P_i*Y;
    R_d=R_i-R_0;
end

T_user=T_set;

diff_count_SSMP=setdiff(T_user,TAC_right);


detect_bit_SSMP=[];
L_index=[];
S_loc_all=[];

  for ii=1:length(T_user)
    index_ii=T_user(ii);
    L_index=[L_index, index_ii];
  end
  
[V_f,l_f]=sort(L_index,'ascend');         %V_f是sort之后的结果，l_f是对应元素标号的顺序
TAC_final=V_f;
S_loc=pinv(H(:,TAC_final))*Y;             %S_loc是sort之后对应的符号向量
S_loc_reshape=reshape(S_loc,[],1);
detect_bit_SSMP=demodulation(S_loc_reshape,log2(M(1)));

end