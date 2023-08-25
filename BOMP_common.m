%2020.11.17 Block-OMP

function [TAC_i,S_i]=BOMP_common(H,M,Y,nTx_AS,nTx,nRx,N0,timeslot) 
%M ���ƽ���

T_i=[];           %��ʼ��
TAC_i=[];
V_s=[];
R_i=Y;

D=kron(H,eye(timeslot));   %������
p=reshape(Y.',[],1);

for i=1:1:nTx_AS
    R_re=H'*R_i;
    R_re=R_re.';
    Dist_R=[];

    for j=1:1:nTx
        B=norm(R_re(:,j),'fro')^2;
        Dist_R=[Dist_R,B];  
    end
    [V,l]=sort(Dist_R,'descend');
    diff_i=setdiff(l,TAC_i,'stable');
    l_i=diff_i(1);
    V_s=[V_s,V(1)];
    TAC_i=[TAC_i l_i];
    [TAC_i,v]=sort(TAC_i,'ascend');    %��������������

%LS Estimation
    W=pinv(H(:,TAC_i));
    H_i=H(:,TAC_i);
    S_i=W*R_i;      %��Ӧ�û������µ��û�����
    S_i_de=reshape(S_i.',[],1);
    D_mod=rx_demodulation(S_i_de,M(1));
    S_i=tx_modulation(D_mod,M(1));
    S_i=reshape(S_i,timeslot,[]);
    S_i=S_i.';
    R_i=R_i-H_i*S_i;
end

