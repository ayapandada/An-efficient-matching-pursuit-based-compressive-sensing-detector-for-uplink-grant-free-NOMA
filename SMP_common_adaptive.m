function  [detect_bit_SSMP,N_error,S_loc,TAC_final,S_Set]=SMP_common_adaptive(H,Y,nRx,nTx,N0,M,nTx_AS,TAC_right,vth,timeslot)
r=nRx-1; %为什么要减1？ 后面用在循环了，感觉这一步没有什么用


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
T_test=[];

R_0=Y;

A=[1:1:nTx];    
D_all=[];

D=1000;
dist=zeros(1,nTx);
Dist=[];                                                             %correction mechinism
 for i_iter=1:nTx_AS    
    T_set_left=setdiff(A,T_set,'stable');
    
    S_left=nTx-length(T_set); %剩余支撑集的大小；
    dist(T_set)=0;
    
    for i1=1:S_left
        k=T_set_left(i1);          
        H_k_i=P_i*H(:,k);
        dist(k)=norm(U_y'*H_k_i)^2/(norm(H_k_i)^2);
    end
    [v,k_hat]=max(dist);           
    [v_test,k_test]=sort(dist);
    
    T_set=[T_set,k_hat];


    B=H(:,T_set)*pinv(H(:,T_set));
    P_i=eye(nRx)-B;
    R_i=P_i*Y;
    diff_R=R_0-R_i;
    D=norm(diff_R,'fro');                                               

    Dist=[Dist,D];
    R_0=R_i;

    if D>=vth*timeslot*N0
     T_test=[T_test,k_hat];
    else
        break;
    end

       
end

T_user=T_test;
L_index=[];
  for ii=1:length(T_user)
    index_ii=T_user(ii);
    L_index=[L_index, index_ii];
  end

  
[V_f,l_f]=sort(L_index,'ascend'); %V_f是sort之后的结果，l_f是对应元素标号的顺序
TAC_final=V_f;

S_loc=pinv(H(:,TAC_final))*Y;             %S_loc是sort之后对应的符号向量
S_loc_reshape=reshape(S_loc,[],1);
detect_bit_SSMP=demodulation(S_loc_reshape,log2(M(1)));
