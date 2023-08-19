clc
clear
G = linspace(0,1000,1000);
H = linspace(0,1000,1000);
FG = linspace(0,1000,1000);
R = linspace(0,1000,1000);
aa = linspace(1,1000,1000);
G(1)=20000;
H(1)=20;
FG(1)=0;
R(1)=0;
ch = 0.2;
SA = 2;
fangcha = 0;
F = 2 + randn(100,1)*2^0.5;
%F_o,SA_o ,pop_1 = shehui(F,FG,G,H,100,SA);
W = 5 + randn(length(G),1)*2^0.5;
for j = 1:100
    a_1 = rand(1);
    if 0.5*(1+tanh(F(j)))>a_1
        X(1,j) = 1;
    else
        X(1,j) = 0;
    end
end
%平均决策
    for j = 1:100
        if F(j) ~= 0
            X(2,j) = exp(tanh(F(j)))/(exp(tanh(F(j)))-1)-1/tanh(F(j));
        else
            X(2,j) = 1/2;
        end
        X(2,j) = (X(2,j)-0.4684)/(0.5816-0.4684);
        X(2,j) = 1*0.5*(1+tanh(F(j)));
    end
    
Gmax=50000; 
record = ones(100,1000);
record_1 = ones(1,1000);
%i = 2;
for i =2:1000   
%% 风险
    if W(i-1)+ch*H(i-1)>H(i-1)
        f = 1-exp(-(W(i-1)+ch*H(i-1)-H(i-1)));
        rt = 1;%记录该次洪水是否造成了农业系统损失
    else
        f = 0;
        rt = 0;
    end
 %% 洪水损失
    if W(i-1)+ch*H(i-1)>H(i-1)
        FG(i-1) = f*(G(i-1)/Gmax)*G(i-1);
    else
        FG(i-1) = 0;
    end    
%% 抬高高度
    if f > 0
        R(i-1) = 1.1*(W(i-1)+ch*H(i-1)-H(i-1));
        if FG(i-1) > 0.5*R(i-1)*(G(i-1)^0.5)
            if G(i-1)-FG(i-1) > 0.5*R(i-1)*(G(i-1)^0.5)
                R(i-1) = 1.1*(W(i-1)+ch*H(i-1)-H(i-1));
            else
                R(i-1) = 0;
            end
        else
            R(i-1) = 0;
        end
    else
        R(i-1) = 0;
    end
 %% 圩堤高度
 fg = FG(i-1);
 g = G(i-1);
 h = H(i-1);
 if i<300
     adv = 0;
 elseif and(i>300,i<400)
     adv = 0.01;
 elseif and(i>400,i<500)
     adv = 0.02;
 elseif and(i>500,i<600)
     adv = 0.05;
 elseif and(i>600,i<700)
     adv = 0.07;
 elseif and(i>700,i<800)
     adv = 0.07;
 elseif and(i>800,i<900)
     adv = 0.09;
 elseif and(i>900,i<1000)
     adv = 0.1;
 end
  if i<300
     adv = 0;
  elseif and(i>300,i<500)
      adv = 2;
  end
    [F,SA ,pop_1,X_o,fangcha_out] = shehui(F,X,fg,g,h,100,SA,fangcha,adv);
    fangcha = fangcha_out;
    F = F';
    record_1(i) =SA;
    record(:,i) = F;
    X = X_o;
    Rm = 0.002*H(i-1)/(1+exp(100-pop_1));
    H(i) = H(i-1) + rt*R(i-1)-rt*0.1*H(i-1) - 0.001*H(i-1) + Rm ;
 %% 农业规模
    G(i) = G(i-1)+ 0.02*G(i-1)*(1-G(i-1)/Gmax)-rt*(FG(i-1)+0.5*R(i-1)*(G(i-1)^0.5));
end
plot(record')
%plot(record_1)
%plot(W)
%hold on
%plot(H)
