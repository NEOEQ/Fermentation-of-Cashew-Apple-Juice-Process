
tic()

//'C' is the vector of concentrations, in which: 'C(1)' viable cells concentration; 'C(2)' dead cells concentration; 'C(3)' substrate concentration and 'C(4)' ethanol concentration, all in g/L
global C

function y=fermenter(V)
    
    //'f' is the derivative of 'C' at a given time
    global f
    global C
    
    // 'to' is the initial time and 'Co' is 'C' at time 'to'
    //'V' is the vector of decision variables, in which: 'V(1)' batch reactor temperature in oC; 'V(2)' initial concentration of substrate in g/L and 'V(3)' operating time in s
    to=0
    Co=[5;0;V(2);0] 
    
    //function containing the kinetic model
    function [f]=fun(t,C,V)
 
        global f
        Ks=23.587; Ki=103.502; Pmax=201.443; ms=0.272; kd=0.054;betamix=2.201; n=0.365; a=2.6; b=-0.282

        mumax=2.305*exp(-61.786/V(1))+-15.326*exp(-152.713/V(1))
        Yxs=-0.187*exp(-83.083/V(1))+0.05*exp(-3.288/V(1))
        Ypx=41.086*exp(-52.601/V(1))+43.330*exp(-52.015/V(1))
        Imix=150
       
        f(1)= (mumax*C(3)/(Ks+C(3)+(C(3))^2/Ki)*(1-C(4)/Pmax))*C(1)-kd*C(1)
        f(2)=kd*C(1)
        f(3)= -((betamix*(1-a*exp(b*C(1)))*C(3)*Imix^n)*((mumax*C(3)/(Ks+C(3)+(C(3))^2/Ki)*(1-C(4)/Pmax))/Yxs +ms))/((betamix*(1-a*exp(b*C(1)))*C(3)*Imix^n)+(((mumax*C(3)/(Ks+C(3)+(C(3))^2/Ki)*(1-C(4)/Pmax))/Yxs +ms)))*C(1)
        f(4)= Ypx*(mumax*C(3)/(Ks+C(3)+(C(3))^2/Ki)*(1-C(4)/Pmax))*C(1)

    endfunction 

    dlist=list(fun,V)
    t = V(3)
    C=ode(Co,to,t,dlist)
    
    //'t_ss' is the time to reach steady state 
    if f(4)<=0.001 then
        f(4)=50
        t_ss=5
        while f(4)>0.001
            t_ss=t_ss+0.1
            C=ode(Co,to,t_ss,dlist)
        end
    else
        t_ss=13
    end
    
    //'w1', 'w2' and 'w3' are the weighting factors --> Weighted Sum Method
    w1=0.3; w2=0.6; w3=0.1
    
    //'dev1' and 'dev2' are penalty terms
    dev1= 32- C(4)    
    dev2=t-t_ss
    
    //objective function
    y=-(w1*(C(4)/t)/9.3 +w2/0.511.*(C(4)./(V(2))) + w3/0.511*(C(4)/(V(2)-C(3))))+ 0.5*abs(max(0,dev1))+ 0.5*abs(max(0,dev2))
endfunction       

////////////////////////////
//PSO algorithm

xmin=[26,70,1]
xmax=[42,170,12]
H = abs(xmax-xmin)
N=30
c1=2
c2=2
D=3

//Inicialization
for i=1:N
    for d=1:D
       x(i,d)=xmin(d)+(xmax(d)- xmin(d))*rand(1,1,'uniform')
       v(i,d)=0
    end
    pk(i,:)=x(i,:)
    P(i)=fermenter(pk(i,:))
    
    if i==1 then
        g=pk(1,:)
        G=P(i)
    end

    if P(i)<G then
       g=pk(i,:)
       G=P(i)
    end
end

subplot(221)
param3d(x(:,1),x(:,2),x(:,3),85,80,"T(°C)@So(g/L)@t(h)")
title("Initial Estimative")

// Set marks to * style
get("hdl").mark_style = 10;

// Turn lines off so we just have points
gce().line_mode="off"
gce().mark_mode="on";
gca().hidden_axis_color = 1

//Iteration

tmax=50

for tatual=1:tmax
    disp(tatual)
    for i=1:N
        R=rand(2,1,'uniform')
        wmax=0.9
        wmin=0.4
        w= wmax-(wmax-wmin)*tatual/tmax
        v(i,:)=w*v(i,:)+ c1*R(1)*(pk(i,:)-x(i,:))+c2*R(2)*(g-x(i,:))
        x(i,:)=x(i,:)+v(i,:)
        for d=1:D
              if x(i,d)<xmin(d)
                  x(i,d)=xmin(d)
                  v(i,d)=0
              end
              if x(i,d)>xmax(d)
                  x(i,d)=xmax(d)
                  v(i,d)=0
              end
        end
        X(i)=fermenter(x(i,:))
        if X(i)<P(i)
            pk(i,:)=x(i,:)
            P(i)=X(i)
            if P(i)< G
                g=pk(i,:)
                G=P(i)
            end
        end
    end
    if tatual==20 then
        subplot(222)
        param3d(x(:,1),x(:,2),x(:,3),85,80,"T(°C)@So(g/L)@t(h)")
        title("Iteration 20")

       // Set marks to * style
        get("hdl").mark_style = 10;

        // Turn lines off so we just have points
        gce().line_mode="off"
        gce().mark_mode="on";
        gca().hidden_axis_color = 1
        
    end
    
    if tatual==40 then
        subplot(223)
        param3d(x(:,1),x(:,2),x(:,3),85,80,"T(°C)@So(g/L)@t(h)")
        title("Iteration 40")

       // Set marks to * style
        get("hdl").mark_style = 10;

        // Turn lines off so we just have points
        gce().line_mode="off"
        gce().mark_mode="on";
        gca().hidden_axis_color = 1
    end
    
    if tatual==50 then
        subplot(224)
        param3d(x(:,1),x(:,2),x(:,3),85,80,"T(°C)@So(g/L)@t(h)")
        title("Iteration 50")

       // Set marks to * style
        get("hdl").mark_style = 10;

        // Turn lines off so we just have points
        gce().line_mode="off"
        gce().mark_mode="on";
        gca().hidden_axis_color = 1
    end
    
end
////////////////////////////////

y=fermenter(g)
x_opt=g

Productivity=(C(4)/x_opt(3))/9.3
Yield1=(1/0.511)*(C(4)/x_opt(2))    
Yield2=(C(4)/(x_opt(2)-C(3)))/0.511

computacional_time=toc()
