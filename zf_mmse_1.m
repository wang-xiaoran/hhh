clear all
clc
% NT=[8 16 32];
% K0=[8 16 32];
nt = 8 ;  %����������
K = 8 ;  %�û��� ����Ϊ������
SNR_V_db = -10:5:30;

SNR_V  = 10.^(-SNR_V_db/10) ;
Error_bits = zeros(2,length(SNR_V)) ; %�����������ʼֵ
yzf=zeros(K,1);
ymmse = zeros(K,1);
loop = 100000; % ѭ������
for k=1:loop  % ѭ������
    bit_source =   randi([0,1],nt,1);% ��Դ; % ���������������
    Xmod = qammod(bit_source,2); %  BPSK ����Ҫ��һ��
    for i=1:length(SNR_V)                
        Noise = 0.707*(randn(K,1)+1j*randn(K,1))*SNR_V(i)^0.5;% ����       
        H=(randn(nt,K)+1j*randn(nt,K))/sqrt(2); % �ŵ����� 
        %%  zf  %%%%%%%%%%
        G = H*H';
        P = H'*pinv(G);
        belta_zf=sqrt(1/trace(P*P'));
        y0zf =belta_zf*H*P*Xmod + Noise;
        for j1= 1:K
            if real(y0zf(j1,1))>=0
                yzf(j1,1) = 1;
            else 
                yzf(j1,1) = 0;
            end
        end
%         yzf=qamdemod(y0zf,2) ;        
        Error_bits(1,i)=Error_bits(1,i)+sum(yzf~=bit_source);       
        %%  mmse  %%%%%%%%%%%%       
        %         W=conj(H)*pinv(H*conj(H)+var(n)*eye(nt)) ;
        %         beta2=sqrt(K/trace(W*conj(W))) ;
        %         y0mmse=beta2*Xmod ;
        %         y1mmse=awgn(y0mmse,n) ;      
        W=H'*inv(H*H'+SNR_V(i)*eye(K));
        belta_mmse = sqrt(1/trace(W*W')) ;
        ylmmse = belta_mmse*H*W*Xmod + Noise;
         for j1= 1:K
            if real(ylmmse(j1,1))>=0
                ymmse(j1,1) = 1;
            else 
                ymmse(j1,1) = 0;
            end
        end
%         ymmse = qamdemod( ylmmse, 2) ;
        Error_bits(2,i) = Error_bits(2,i) + sum( ymmse ~= bit_source ) ;        
    end    
end
BER = Error_bits/(loop*nt);
BER_zf = BER(1,:);
BER_mmse = BER(2,:);
figure(1);
semilogy(SNR_V_db,BER_zf,'-b*');
hold on
% semilogy(SNR_V_db,BER_mmse,'-ro');
% legend('ZF','MMSE');
xlabel('SNR');
ylabel('BER');