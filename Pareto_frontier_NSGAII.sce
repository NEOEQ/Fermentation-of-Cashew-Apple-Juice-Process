tic()

function g=fermenter(V)
    
    //'C' is the vector of concentrations, in which: 'C(1)' viable cells concentration; 'C(2)' dead cells concentration; 'C(3)' substrate concentration and 'C(4)' ethanol concentration, all in g/L
    //'f' is the derivative of 'C' at a given time
    global f
    
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

    t = V(3)
    dlist=list(fun,V)
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
    
    //'dev1' and 'dev2' are penalty terms
    dev1= 32- C(4)
    dev2=t-t_ss
    
    //objetive functions
    g(1)=-(C(4)/t)/9.3 +0.5*(max(0,dev1))+ 0.5*abs(max(0,dev2))//  Productivity
    g(2)=-(1/0.511)*(C(4)/V(2)) +0.5*(max(0,dev1))+ 0.5*abs(max(0,dev2))//  Yield1
    g(3)=-(C(4)/(V(2)-C(3)))/0.511 +0.5*(max(0,dev1))+ 0.5*abs(max(0,dev2))//  Yield2


endfunction

//NSGA-II’s specifications
PopSize     = 506;
Proba_cross = 0.5;
Proba_mut   = 0.4;
NbGen       = 100;
NbCouples   = 110;
Log         = %T;
nb_disp     = 10; // Nb point to display from the optimal population
pressure    = 0.1;

ga_params = init_param();
ga_params = add_param(ga_params,'dimension',3);
ga_params = add_param(ga_params,'minbound',[26,70,0]);
ga_params = add_param(ga_params,'maxbound',[42,170,12]);

//NSGA-II algorithm
[pop_opt, fobj_pop_opt, pop_init, fobj_pop_init] = optim_nsga2(fermenter, PopSize,NbGen, Proba_mut, Proba_cross, Log, ga_params)

computacional_time=toc()

for b=1:PopSize
    pop_init_1(b,:)=pop_init(b)
    pop_opt_1(b,:)=pop_opt(b)
end


//Plotting Pareto frontier and decision variables
scf()
scatter3(-fobj_pop_opt(2:$,3),-fobj_pop_opt(2:$,1),-fobj_pop_opt(2:$,2),'fill')
gca().rotation_angles = [60, 45];
xlabel('$\eta_{P/\Delta S}$')
ylabel('$\frac{Q_{P}}{Q_{P}_{max}}$')
zlabel('$\eta_{P/S_{o}}$')
gca().hidden_axis_color = 1

scf()
scatter3(pop_opt_1(2:$,1),pop_opt_1(2:$,2),pop_opt_1(2:$,3),'fill')
gca().rotation_angles = [60, 45];
xlabel('$T(^{o}C)$')
ylabel('$S_{o}(g/L)$')
zlabel('$t_{F}(h)$')
gca().hidden_axis_color = 1
