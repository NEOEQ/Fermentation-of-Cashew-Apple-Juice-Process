
//function containing the kinetic model
function [f]=fun(t,C,T)
        
        Ks=23.587; Ki=103.502; Pmax=201.443; ms=0.272; betamix=2.201; n=0.365; a=2.6; b=-0.282;kd=0.054;
        mumax=2.305*exp(-61.786/T)-15.326*exp(-152.713/T)
        Yxs=-0.187*exp(-83.083/T)+0.05*exp(-3.288/T)
        Ypx=+41.086*exp(-52.601/T)+43.330*exp(-52.015/T)
        Imix=150
        
        f(1)= (mumax*C(3)/(Ks+C(3)+(C(3))^2/Ki)*(1-C(4)/Pmax))*C(1)-kd*C(1)
        f(2)=kd*C(1)
        f(3)= -((betamix*(1-a*exp(b*C(1)))*C(3)*Imix^n)*((mumax*C(3)/(Ks+C(3)+(C(3))^2/Ki)*(1-C(4)/Pmax))/Yxs +ms))/((betamix*(1-a*exp(b*C(1)))*C(3)*Imix^n)+(((mumax*C(3)/(Ks+C(3)+(C(3))^2/Ki)*(1-C(4)/Pmax))/Yxs +ms)))*C(1)
        f(4)= Ypx*(mumax*C(3)/(Ks+C(3)+(C(3))^2/Ki)*(1-C(4)/Pmax))*C(1)

endfunction 

// batch reactor temperatures in oC
T_a=30 //(a)
T_b=30  //(b)
T_c=34 //(c)
T_d=42 //(d)

to=0 //initial time

//vectors of concentrations at initial time: Co(1)= viable cells concentration; Co(2)= dead cells concentration; Co(3)= substrate concentration and Co(4)= ethanol concentration, all in g/L
Co_a=[4.21;0;168.83;0] //(a)
Co_b=[3.70;0;72.54;0] //(b)
Co_c=[5.33;0;106.85;0] //(c)
Co_d=[5;0;101.65;0] //(d)

//vector of time in h
t=0:0.1:12 

dlist_a=list(fun,T_a)
dlist_b=list(fun,T_b)
dlist_c=list(fun,T_c)
dlist_d=list(fun,T_d)

////vectors of concentrations at time 't' in g/L
C_a=ode(Co_a,to,t,dlist_a) //(a)
C_b=ode(Co_b,to,t,dlist_b) //(b)
C_c=ode(Co_c,to,t,dlist_c) //(c)
C_d=ode(Co_d,to,t,dlist_d) //(d)

//experimental data in g/L

//(a)
X_exp_a=[4.21,4.88,5.22,6.65,7.49,7.91,8.24]
S_exp_a=[168.83,154.83,116.97,100.50,62.08,45.05,24.42]
P_exp_a=[0.00,6.14,14.38,31.21,40.04,56.45,63.18]

//(b)
P_exp_b=[0.00,11.05,24.35,30.52,30.97,30.74,31.19]
S_exp_b=[72.54,46.20,17.84,3.87,0.00,0.00,0.00]    
X_exp_b=[3.70,4.66,5.44,5.78,6.12,5.78,5.83]

//(c)
P_exp_c=[0, 20.422535,37.840374,50.985916,51.643192,52.957745]
S_exp_c=[106.85446,73.333336,33.568073,7.6056337,0,0]
X_exp_c=[5.3051643, 6.619718, 7.276995, 7.6056337, 7.6056337, 7.6056337]

//(d)
P_exp_d=[0.00,17.851852,33.482105,43.925926,44.5,45.68]
S_exp_d=[101.65,76.19,38.96,8,0.00,0.00]
X_exp_d=[5.00,6.75,7.18,7.72,7.53,7.66]


t_exp_a=[ 0, 2, 4, 6, 8, 10,12] //(a)
t_exp_b=[ 0, 2, 4, 6, 8, 10,12] //(b)
t_exp_c=[ 0, 2, 4, 6, 8, 10] //(c)
t_exp_d=[ 0, 2, 4, 6, 8, 10] //(d)

//plotting experimental data and model
subplot(221)
X_a= C_a(1,:)+C_a(2,:)
plot(t,X_a,'b',t,C_a(3,:),'g',t,C_a(4,:),'r',t_exp_a,X_exp_a,'bo',t_exp_a,S_exp_a,'gx',t_exp_a,P_exp_a,'rs')  
xtitle('(a)','','Concentration (g/L)')
   
subplot(222)
X_b= C_b(1,:)+C_b(2,:)
plot(t,X_b,'b',t,C_b(3,:),'g',t,C_b(4,:),'r',t_exp_b,X_exp_b,'bo',t_exp_b,S_exp_b,'gx',t_exp_b,P_exp_b,'rs') 
xtitle('(b)','','')

subplot(223)
X_c= C_c(1,:)+C_c(2,:)
plot(t,X_c,'b',t,C_c(3,:),'g',t,C_c(4,:),'r',t_exp_c,X_exp_c,'bo',t_exp_c,S_exp_c,'gx',t_exp_c,P_exp_c,'rs')   
xtitle('(c)','Time (h)','Concentration (g/L)')

subplot(224)
X_d= C_d(1,:)+C_d(2,:)
plot(t,X_d,'b',t,C_d(3,:),'g',t,C_d(4,:),'r',t_exp_d,X_exp_d,'bo',t_exp_d,S_exp_d,'gx',t_exp_d,P_exp_d,'rs')  
legend(['Biomass','Substrate','Ethanol'],1)  
xtitle('(d)','Time (h)','')

    
       
