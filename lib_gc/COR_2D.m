function [sync,async] = COR_2D(T,data)
%�ȼ���Ƿǵȼ����������������Ա�����Ҳ�����ڵȼ���Ŷ��Ĺ���
[n,m] = size(data);
noda  = zeros(n,n);
%% �¶���չ
t0 = 2*T(1)-T(2);
tend = 2*T(end)-T(end-1);
T_expand = [t0;T(:);tend];
%% noda�������
for i=1:n
    for j=1:n
        if i~=j
            noda(i,j) = (T_expand(j+2)-T_expand(j))/(2*pi)/(T_expand(j+1)-T_expand(i+1));  %noda��������㹫ʽ
        end
    end
end
%% ��ȡ��̬����
sum_data = zeros(1,m);
deltaT_sum = 0;
for i = 1:n
    sum_data = sum_data+(T_expand(i+2)-T_expand(i))*data(i,:);
    deltaT_sum = deltaT_sum+(T_expand(i+2)-T_expand(i));
end
mean_data = sum_data/deltaT_sum;
data = data-mean_data(ones(n,1),:);
%% ͬ���׺��첽�׼���
weightT = T_expand(3:end) - T_expand(1:end-2);
weightT = diag(weightT);
sync = zeros(m,m);
async = zeros(m,m);
for i = 1:m
    for j = 1:m
        sync(i,j) = data(:,i)'*weightT*data(:,j);
        async(i,j) = data(:,i)'*weightT*noda*data(:,j);
    end
end
sync = sync/(2*(T(end)-T(1)));
async = async/(2*(T(end)-T(1)));
end
