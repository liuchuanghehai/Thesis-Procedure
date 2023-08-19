function [F_o,SA_o ,pop_1,record_3,fangcha_out] = shehui(F,X,FG,G,H,L,SA,fangcha_in,adv_in)
%% 基本参数
sizepop = 100;
ger = 1;
%r_1 = 0+3*rand(sizepop,1);
%r_1 = 0 + randn(sizepop,1)*2^0.5;
r_1 = F;
for iiii = 1:length(r_1)
    if r_1(iiii,1)>3
        r_1(iiii,1) = 3;
    elseif r_1(iiii,1)<-3
        r_1(iiii,1) = -3;
    end
end
x_s = linspace(1,10,10);
y_s = linspace(1,10,10); 
%% 初始化输入
site = ones(2,100);
k_r = 0;
for k_x = 1:10   
    for k_y = 1:10
        site(1,k_r*10+k_y) = x_s(k_x);%x坐标
        site(2,k_r*10+k_y) = y_s(k_y);%y坐标
    end
   
    k_r = k_r + 1;
end
for j = 1:sizepop
    pop(1,j) = site(1,j);%x坐标
    pop(2,j) = site(2,j);%y坐标
    pop(3,j) = r_1(j,:)';%意识t
    if pop(3,j)>3
        pop(3,j)=3;
    elseif pop(3,j)<-3
        pop(3,j)=-3;
    end
    pop(4,j) = 0;%行为t
    pop(5,j) = 0;%农业规模t
    pop(6,j) = 0;%均值行为t
    pop(7,j) = 0;%均值规模t
    pop(8,j) = 0;%意识t+1
    pop(9,j) = 0;%行为t+1
    pop(10,j) = 0;%农业规模
    pop(11,j) = 0;%平均行为
    pop(12,j) = 0;%平均规模
end
%% 循环
record_1 = ones(ger,1);
record_2 = ones(ger,sizepop);
record_3 = ones(2,sizepop);
pop(9,:) = X(1,:);
pop(11,:) = X(2,:);
k = 1;
fangcha = fangcha_in;
while k <=ger
if FG~=0
    rrr = 1;
else
    rrr = 0;
end
%% 收益

    for j = 1:sizepop
        a_2 = rand(1);
        if a_2 >0
            cf = 0.5;
        else 
            cf = 0;
        end
        %收益
        if pop(9,j) == 1
            pop(10,j) = (G- FG)/L-0.5*0.001*H*(G^0.5)/L-0.3*G/L;
            pop(12,j) = (G- FG)/L-0.5*0.001*H*(G^0.5)/L*pop(11,j)-0.3*G*pop(11,j)/L-rrr*cf*G*(1-pop(11,j))/L;
        elseif pop(9,j) == 0
            pop(10,j) = (G- FG)/L-rrr*cf*G/L;
            pop(12,j) = (G- FG)/L-0.5*0.001*H*(G^0.5)/L*pop(11,j)-0.3*G*pop(11,j)/L-rrr*cf*G*(1-pop(11,j))/L;
        end
    end
        %平均收益            
%% 意识
    heb_aiout = 0;
    heb_a_heb = 0;
    for j = 1:sizepop
        %% 平均场
        ex = pop(11,j);
        x = pop(9,j);
        paix = pop(10,j);
        paiex = pop(12,j);
        j_x = pop(1,j);
        j_y = pop(2,j);
        heb = 0;
        heb_a = 0;
        heba = 0;
        if max((paix-paiex)/paiex,0)>0
            kk2 = 1;
        else
            kk2 = 0;
        end
        %kk2 = max((paix-paiex)/paiex,0);        
        for i = 1:sizepop
            i_x = pop(1,i);
            i_y = pop(2,i);
            d = ((i_x-j_x)^2+(i_y-j_y)^2)^0.5;
            kk1 = max(1-abs(ex-pop(9,i))/1,0);
            fai = 100*kk1*kk2*((pop(10,i)-paiex)/paiex);              
            if d == 0
                b = 0;
            else
                b = fai/d;
            end
            if b < 0.0
                b = 0;
            end
            heb = heb + b;
            if pop(9,i) == 1
                heba = heba + b*(3);
            elseif pop(9,i) == 0
                heba = heba + b*(-3);
            end

            if max((pop(10,i)-pop(12,i))/pop(12,i),0)>0
                kk2_a = 1;
            else
                kk2_a = 0;
            end
            %kk2_a = max((pop(10,i)-pop(12,i))/pop(12,i),0);
            kk1_a = max(1-abs(pop(11,i)-x)/1,0);
            fai_a = 100*kk1_a*kk2_a*((paix-pop(12,i))/pop(12,i));               
            if d == 0
                b_a = 0;
            else
                b_a = fai_a/d;
            end
            if b_a <0.0
                b_a = 0;
            end
            heb_a = heb_a + b_a;
        end%对意识内部交流结束
            if heb~=0
                Aba = heba/heb;
                da1 = 0.1*(Aba-pop(3,j));
            else
                da1 = 0;
            end
            %逆作用，推球社会意识
            heb_aiout = heb_aiout + heb_a*pop(3,j);
            heb_a_heb = heb_a_heb + heb_a;         
        %% 社会意识
        da2 = (1/(1+exp(-abs(SA-pop(3,j))))-0.5)*(SA-pop(3,j));
        %% 外势能场
        adv = adv_in;
        if fangcha~=0
            dou = (1-heb_a/fangcha)*(adv);
        else
            dou = 0;
        end
        da3 = dou*(3-pop(3,j));
        if da3 ~=0
            da3;
        end      
        pop(8,j) = pop(3,j)+ 0.3*da1 + 0.01*da2 + 0.2*da3;
        if pop(8,j)>3
            pop(8,j) = 3;
        elseif pop(8,j)<-3
            pop(8,j) = -3;
        end 
        record_2(k,j) = pop(8,j);       
    end %对意识演化
    fangcha_out = heb_a_heb;
    if heb_a_heb~=0
        Ainba = heb_aiout/heb_a_heb;
        dsa = 0.2*(1-1/(1+exp(-abs(Ainba-SA))))*(Ainba-SA);
    else
        dsa = 0;
    end     
    SA = SA + dsa;
    if SA >3
        SA = 3;
    elseif SA<-3
        SA = -3;
    end    
%% 决策    
%确定的决策
pop_1 = 0;
for j = 1:sizepop
    a_1 = rand(1);
    if 0.5*(1+tanh(pop(3,j)))>a_1
        pop(9,j) = 1;
        pop_1 = pop_1+1;
    else
        pop(9,j) = 0;
    end
end
%平均决策
    for j = 1:sizepop
        if pop(3,j) ~= 0
            pop(11,j) = exp(tanh(pop(3,j)))/(exp(tanh(pop(3,j)))-1)-1/tanh(pop(3,j));
        else
            pop(11,j) = 1/2;
        end
        pop(11,j) = (pop(11,j)-0.4684)/(0.5816-0.4684);
    end
    pop(11,j) = 1*0.5*(1+tanh(pop(8,j)));
record_1(k) = pop_1(1);
k = k+1;
end
record_3(1,:) = pop(9,j);
record_3(2,:) = pop(11,j);
F_o = record_2;
SA_o = SA;
end
%plot(record_1);title('执行数量迭代')
%hold on
%plot(record_2)