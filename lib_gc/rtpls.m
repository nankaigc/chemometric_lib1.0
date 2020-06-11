%% ����ѡ�񷽷�֮RT���������������
function [result]=rtpls(cal, caltar, nLV,n_mc, selectLV)  %����cal��У��������ģ����caltar��У������Ӧ��Ũ�ȣ�runtimes�ǵ���������һ��ѡ��500��
%factor��ѡ�����������һ������һ���淨LOOV��ã����ǿ�����

method='center';
n_fold=5;
order=0;

tic

[num,col]=size(cal);
cal=cal-mean(cal);
caltar=caltar-mean(caltar);


B = simpls(cal, caltar, nLV);
coef0 = B(:,end);

for i=1:n_mc
    rlist = randperm(num);
    y2 = caltar(rlist);
    B = simpls(cal, y2, nLV);
    coefs(i,:) = B(:,end)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a,b] = size(coefs);  P=[];
for i = 1:length(coef0)
    noisev = coef0(i);
    row = find(abs(coefs(:,i))>=abs(noisev));
    k = length(row);
    per = k/a;
    P = [P,per];
end

[~,index] = sort(P);  %��ʱpp���ǰ��ս������еģ�Ȼ��index���Ƕ�Ӧ�Ĳ����㣬wvleng�еĲ�����

%+++  Cross-Validation to choose an optimal subset;
RMSECV=zeros(1,col);
Q2_max=zeros(1,col);
LV=zeros(1,col);
for i=1:col
    vsel=index(1:i);
    
    CV=plscv(cal(:,vsel),caltar,nLV,n_fold,method,0,order);
    if selectLV == 0
        RMSECV(i)=CV.RMSECV_min;
        Q2_max(i)=CV.Q2_max;
        LV(i)=CV.optLV;
    elseif selectLV==1
        RMSECV(i)=CV.RMSECV_min_1SD;
        Q2_max(i)=CV.Q2_max_1SD;
        LV(i)=CV.optLV_1SD;
    end
    if mod(i,100)==0
        fprintf('The %d/%dth variable validation finished.\n',i,col);
    end
end
[RMSECV_min,indexOPT]=min(RMSECV);
Q2_max=max(Q2_max);


%+++ save results;
time=toc;
%+++ output
result.P_value=P;
result.time=time;
result.RMSECV=RMSECV;
result.RMSECV_min=RMSECV_min;
result.Q2_max=Q2_max;
result.iterOPT=indexOPT;
result.optLV=LV(indexOPT);
result.vsel=index(1:indexOPT);
end


%�����ǿ�P�Ĵ�С������P���ǰ���˳������ģ���Ҫ��P���н�������
% [pp,index] = sort(P,'descend');  %��ʱpp���ǰ��ս������еģ�Ȼ��index���Ƕ�Ӧ�Ĳ����㣬wvleng�еĲ�����
% ind = index(1:220);   %һ��ѡ���ܲ�����10%��20%����Ҫ����20%��һ��2205�������㣬��������ѡ����220��
% sp = data.sp(:,ind);  %���������ѡ��Ĳ�����