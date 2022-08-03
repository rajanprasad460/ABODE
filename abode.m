function abode(sys,y)
%---------Finding poles and zeros --------------
if nargin == 1
    k=1;
elseif nargin == 2
    k=y;
else
    error('arguments must be either ''abode(sys,k)'' -- or ''abode(sys)''');
end


clc;
%--------Finding hidden constant---------------------
[num,den]=tfdata(sys,'v');
k2=1/den(1);
k3=find(num>0);
k3=num(k3(1))*k;
syms s;
tit=poly2sym(num,s)/poly2sym(den,s);
tit=simplify(tit);
stit=sprintf('Bode Plot for TF %s',tit);
%------------Ploes and zeros----------------------------
p_loc=pole(sys);
z_loc=zero(sys);

c_fp=sort(abs(p_loc));
c_fz=sort(abs(z_loc));
c_f=sort([abs(p_loc);abs(z_loc)]);

uc_f(1,:)=unique(c_f);
czp=unique(c_f);

n=length(uc_f);
for i=1:n
    uc_f(2,i)=length(find(c_f==uc_f(1,i)));
end

gslope=zeros(1,n);
tslope=0;
dummy_g=zeros(1,n);

% -----Data for slope for gain asymptotes---------
com_gain=zeros(n,n);
for i=1:n
    temp_c_f=uc_f(1,i);
    gd_z=ismember(temp_c_f,c_fz);
    if gd_z==1
        tslope=tslope+20*uc_f(2,i);
        dummy_g(uc_f(1,:)>=temp_c_f)=20*uc_f(2,i);
    else
        tslope=tslope-20*uc_f(2,i);
        dummy_g(uc_f(1,:)>=temp_c_f)=-20*uc_f(2,i);
    end
    gslope(i)=tslope;
    com_gain(i,:)=dummy_g;
    dummy_g=zeros(1,n);
end
com_gain(end+1,:)=gslope;


%--------------finding poles at zero and assigning slopes---------------
op_slope=[];
oz_slope=[];
aop_con=0;
if (ismember(0,c_fp)==1)
    op_slope=uc_f(2,1)*(-20);
    aop_con=uc_f(2,1)*(-90);
    c_fp(1:uc_f(2,1))=[];
    uc_f(:,1)=[];
    gslope(1)=[];
end
if (ismember(0,c_fz)==1)
    oz_slope=uc_f(2,1)*(20);
    aop_con=uc_f(2,1)*(90);
    c_fz(1:uc_f(2,1))=[];
    uc_f(:,1)=[];
    gslope(1)=[];
end
%--------------------Plotting bode and reading axes----------------

[~,~,w]=bode(k*sys);
bode(k*sys)
h=get(findall(get(gcf,'Children'),'String','Magnitude (dB)'),'Parent');
axes(h);
yxlim=axis;
hold on;
grid;

%------------Plotting constant magnitude line 1st 0 c_f---------
if (op_slope~=0)
    ct=20*(log10((k3*k2*prod(c_fz)/prod(c_fp))/(uc_f(1,1))^(length(c_f)-sum(uc_f(2,:)))));
elseif (oz_slope~=0)
    ct=20*(log10((k3*k2*prod(c_fz)/prod(c_fp))*(uc_f(1,1))^(length(c_f)-sum(uc_f(2,:)))));
else
ct=20*(log10(k3*k2*prod(c_fz)/prod(c_fp)));
end

plot([w(1) w(end)],[ct ct],'--r','LineWidth',1);
%---for poles/zeros at zero-----------------------------------
if(op_slope~=0)
    plot([w(1) uc_f(1,1)],[ct+log10(w(1)/uc_f(1,1))*op_slope ct],'LineWidth',1)
end
if(oz_slope~=0)
    plot([w(1) uc_f(1,1)],[ct+log10(w(1)/uc_f(1,1))*oz_slope ct],'LineWidth',2)
end
%---------------------------------------------------------------------------------

p1=1;
ct_t=ct;
n=length(uc_f(1,:));
while (p1<=n)
    if p1==n
        ct_t=ct_t+log10(w(end)/uc_f(1,p1))*gslope(p1);
        plot([uc_f(1,p1) w(end)],[ct ct_t],'LineWidth',2)
    else
        ct_t=ct_t+log10(uc_f(1,p1+1)/uc_f(1,p1))*gslope(p1);
        plot([uc_f(1,p1) uc_f(1,p1+1)],[ct ct_t],'LineWidth',2)
        ct=ct_t;
    end
    p1=p1+1;
end



%--------------Phase Plot-------------------------------
phasec_f=[uc_f(1,:)*0.1 uc_f(1,:)*10];
dummy_a=zeros(size(phasec_f));
phase_a=zeros(size(phasec_f));
com_phase=zeros(n,2*n);
for i=1:n
    temp_c_f=uc_f(1,i);
    gd_z=ismember(temp_c_f,c_fz);
    
    if gd_z==1
        dummy_a(phasec_f>=temp_c_f*0.1 & phasec_f<temp_c_f*10)=45*uc_f(2,i);
    else
        dummy_a(phasec_f>=temp_c_f*0.1 & phasec_f<temp_c_f*10)=-45*uc_f(2,i);
    end
    phase_a=phase_a+dummy_a;
    com_phase(i,:)=dummy_a;
    dummy_a=zeros(size(phasec_f));
end
com_phase(end+1,:)=phase_a;


h=get(findall(get(gcf,'Children'),'String','Phase (deg)'),'Parent');
axes(h);
hold on;
grid on;

plot([w(1) w(end)],[aop_con aop_con],'--r');
apt=aop_con;
for i=1:2*n
    if (i==2*n)
        apt_t=apt+log10(w(end)/phasec_f(i))*phase_a(i);
        plot([phasec_f(i) w(end)],[apt apt_t])
    else
        apt_t=apt+log10(phasec_f(i+1)/phasec_f(i))*phase_a(i);
        plot([phasec_f(i) phasec_f(i+1)],[apt apt_t]);
        apt=apt_t;
    end
end

%-----------Assigning title------------------------
title(stit);
%--------Disp Gain Magnitude table for asymptotes--------
newline;
disp(' Slope for Bode magnitude Asymptotes');

T = array2table(com_gain);

for i=1:length(czp)
    st_nill{i}=sprintf('Start %.1f',czp(i));
    srows{i}=sprintf('w=%.2f ',czp(i));
end
srows{i+1}='Sum of Slopes';
T.Properties.VariableNames=st_nill';
T.Properties.RowNames=srows';
disp(T);

%---------Display phase angle table of asymptotes---------------------
newline;
disp('Slope for Bode phase Asymptotes');


T2 = array2table(com_phase);
st_nill2='';

for i=1:2*n
    if (i<=n)
        st_n='Start';
    else
        st_n='End';
    end
    st_nill2{i}=sprintf('%s %.2f',st_n,phasec_f(1,i));
end
for i=1:n
    psrows{i}=sprintf('w=%.2f ',czp(i));
end
psrows{i+1}='Sum of Slopes';
T2.Properties.VariableNames=st_nill2';
T2.Properties.RowNames=psrows';
disp(T2);
hold off;
end