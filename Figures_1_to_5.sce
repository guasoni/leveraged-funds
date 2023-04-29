
clear
funcprot(0)

sigma=0.16

function f= perf(pim,pip,epsilon,kappa)
mu=kappa*sigma^2


    A=pim/(1-pim)
    B=pip/(1-pip)
    nc=sign(A)/(2*kappa-1)*(abs(B)^(2*kappa-1)-abs(A)^(2*kappa-1))

    mean_plus_ATC=mu*integrate('1/nc*z./(1+z).*((z.^2).^(kappa-1))','z',A,B)
    G=epsilon/(1+B)/(1+(1-epsilon)*B)
    ATC=sigma^2*(2*kappa-1)/2*G*B/(1-(A/B)^(2*kappa-1))
    m=mean_plus_ATC-ATC;
    st_dev=sigma*sqrt(integrate('1/nc*z.^2 ./(1+z).^2.*((z.^2).^(kappa-1))','z',A,B))
    f=[ATC,m,st_dev]
endfunction



function f=integrand(z,A,gam,L,kappa)
    mu=kappa*sigma^2
      
     hz=(gam*sigma^2*L)*(z./(1+z))-gam*sigma^2/2*(z./(1+z)).^2;
     hA=(gam*sigma^2*L)*(A./(1+A))-gam*sigma^2/2*(A./(1+A)).^2;
    f=(hz-hA).*(z/B).^(2*mu/sigma^2-2);
endfunction

function f=WAB(A,B,gam,L,kappa)
    
    f1=integrate('integrand(z,A,gam,L,kappa)','z',A,B);
    f=2./B.^2/sigma^2 * f1;

endfunction

function f=eq1(A,B,gam,epsilon,L,kappa)
   
    
    
    W=WAB(A,B,gam,L,kappa);
    term2=epsilon./(1+B)./(1+(1-epsilon)*B);

    f=W-term2;
endfunction

function f=eq2(A,B,gam,epsilon,L,kappa)
    mu=kappa*sigma^2
     term1a=2/sigma^2./B.^2;
     
     term1b=(gam*sigma^2*L)*(B./(1+B))-gam*sigma^2/2.*(B./(1+B)).^2-((gam*sigma^2*L)*(A./(1+A))-gam*sigma^2/2*(A./(1+A)).^2);
     term1=term1a*term1b;
     
     term2=2*mu/sigma^2./B;
     term2=term2.*WAB(A,B,gam,L,kappa);
     
     term3=(1-epsilon)^2./(1+(1-epsilon)*B).^2-1./(1+B).^2;
    f=term1-term2-term3;
endfunction


function f=freeboundaryproxies(gam,epsilon,L,kappa)
   

   dev1=(3/4/gam*L^2*(L-1)^2*epsilon)^(1/3)
   dev2=(L-kappa/2)/gam*signm(L*(L-1))*(gam/6*abs(L*(L-1)))^(1/3)*epsilon^(2/3)
   
   pim=L-dev1-dev2
   pip=L+dev1-dev2
   
   A=pim/(1-pim)
   B=pip/(1-pip)
   
      
   
f=[A,B];
endfunction



function f=tradingboundaries(A0,B0, gam,epsilon,L,kappa)
        mu=kappa*sigma^2
        //System of equations:
        deff('res=fct(x)',['res(1)=eq1(x(1),x(2),gam,epsilon,L,kappa)';'res(2)=eq2(x(1),x(2),gam,epsilon,L,kappa)']);
        x0=[A0; B0];
      
        AB=fsolve(x0, fct,1.d-4);
        A1=AB(1);
        B1=AB(2);
      
       
pim=A1/(1+A1)
pip=B1/(1+B1)

//Realized Multiplier=int(pi)dt (equals beta for kappa=0)
RM=pim.*pip ./(pip-pim) .*(log(pip ./pim));
//Alpha
alpha=-sigma^2/2*pim .*pip.*(pip-1).^2 ./(pip-pim) ./(1/epsilon-pip);

//Tracking Error proxy in terms of exact alpha
Iofpi2=pim.*pip

TRE=sqrt(sigma^2*Iofpi2-2*L*sigma^2*RM+sigma^2*L^2)



G=epsilon/((1+B1)*(1+(1-epsilon)*B1))
ATC_g=sigma^2/2 *B1 *G*(2*kappa-1)/(1-(A1/B1)^(2*kappa-1))
hA=gam*L*pim-gam/2 *pim^2

termx=sigma^2*hA-sigma^2*(1-2*kappa)/2*gam*pim+ATC_g*(1-gam*(pip-pim)/B1/G)
m_g0x=1/gam/(L-kappa)*termx
m_gx=m_g0x*kappa


termy=2*ATC_g*(kappa-gam*L*(pip-pim)/B1/G)+2*kappa*sigma^2*hA+gam*L*(1-2*kappa)*sigma^2*pim
v_g=1/gam/(L-kappa)*termy
alpha_g=-ATC_g//only true, if mu=0.
tri= perf(pim,pip,epsilon)

m_g0=(sigma^2*hA+ATC_g+gam/2*v_g)/gam/L
m_g=m_g0*kappa



TRE_g=sqrt(v_g-2*L*m_g0+sigma^2*L^2)


mean_g=tri(2)
st_g=tri(3)
f=[gam,pim,pip,ATC_g,alpha_g,TRE_g,mean_g,st_g];
endfunction


//      Levered ETF: Trading Boundaries          //


function f=perf_trading(gam,kappa,L,epsi,produce)

vec1=zeros(1,8);

//Compute proxy for first gamma
epss=1/10000;
o=freeboundaryproxies(gam(1),epss,L,kappa);
A=o(1);
B=o(2);
//iterate until epss equal to actual transaction costs
while epss<epsi
dd=tradingboundaries(A,B,gam(1),epss,L,kappa);
epss=epss*1.2;
A=dd(2)/(1-dd(2));
B=dd(3)/(1-dd(3));
end
//exact performance statistics
dd=tradingboundaries(A,B,gam(1),epsi,L,kappa);
vec1(1,:)=dd;

//Iterate now with decreasing gamma
//This is only relevant, when gamma is a vector

//try 
    
for i=2:length(gam)

A0=vec1(i-1,2)/(1-vec1(i-1,2));
B0=vec1(i-1,3)/(1-vec1(i-1,3));
dd=tradingboundaries(A0,B0,gam(i),epsi,L,kappa)
disp(dd)

vec1(i,:)=dd;
end

//Turn chronological order
vec1=vec1($:-1:1,:);
//catch
//end
//Save performance in File
if produce==1 then

    if L>0 then name="Bull"; end
    if L<0 then name="Bear"; end
    name=name+string(abs(L*10))+"bp"+string(epsi*10000)+"kappa"+string(kappa)
    fprintfMat(name, vec1, "%10.12f");
end

f=vec1
endfunction

//Figure 1 and 2: Vanishing Risk Premium. Multipliers Lambda -3,-2,-1,2,3,4. Large range of gamma.
// Uncomment the following lines so to reproduce files for optimal trading boundaries for Bullish and Bearish funds.
gam=[900000:-100000:100000,99000:-1000:1000,900:-100:100,99.9:-0.1:0.5, 0.499:-0.01:0.09, 0.0899:-0.0001:0.0012,0.001201:-0.00001:0.000051];
perf_trading(gam,0,-3,0.001,1)
gam=[900000:-100000:100000,99000:-1000:1000,900:-100:100,99.9:-0.1:0.5, 0.499:-0.01:0.09, 0.0899:-0.0001:0.0012,0.001201:-0.00001:0.000051];
perf_trading(gam,0,-2,0.001,1)
gam=[900000:-100000:100000,99000:-1000:1000,900:-100:100,99.9:-0.1:0.5, 0.499:-0.01:0.09, 0.0899:-0.0001:0.0012,0.001201:-0.00001: 0.000081];
perf_trading(gam,0,-1,0.001,1)
gam=[900000:-100000:100000,99000:-1000:1000,900:-100:100,99.9:-0.1:0.5, 0.499:-0.01:0.09, 0.0899:-0.0001:0.0012,0.001201:-0.00001: 0.001021000000];
perf_trading(gam,0,2,0.001,1)
gam=[900000:-100000:100000,99000:-1000:1000,900:-100:100,99.9:-0.1:0.5, 0.499:-0.01:0.09, 0.0899:-0.0001:0.0012,0.001201:-0.00001:0.000521000000];
perf_trading(gam,0,3,0.001,1)
gam=[900000:-100000:100000,99000:-1000:1000,900:-100:100,99.9:-0.1:0.5, 0.499:-0.01:0.09, 0.0899:-0.0001:0.0012,0.001201:-0.00001:0.000341000000];
perf_trading(gam,0,4,0.001,1)


//While the above 6 lines are run, trading boundaries and some performance data are saved in textfiles. These are read now:


//Reading Files

dm3ETFbp10=fscanfMat('Bear30bp10kappa0');
dm2ETFbp10=fscanfMat('Bear20bp10kappa0');
dm1ETFbp10=fscanfMat('Bear10bp10kappa0');

dp2ETFbp10=fscanfMat('Bull20bp10kappa0');
dp3ETFbp10=fscanfMat('Bull30bp10kappa0');
dp4ETFbp10=fscanfMat('Bull40bp10kappa0');


//Figure 1
clf

L1=0.01/100:0.01/100:1
L2=L1*0
L3=L2+1
plot(L1,L2,'k',L1,L3,'k',dp4ETFbp10(:,6),dp4ETFbp10(:,3),'r',dp4ETFbp10(:,6),dp4ETFbp10(:,2),'r--',dp3ETFbp10(:,6),dp3ETFbp10(:,3),'b',dp3ETFbp10(:,6),dp3ETFbp10(:,2),'b--',dp2ETFbp10(:,6),dp2ETFbp10(:,3),'g',dp2ETFbp10(:,6),dp2ETFbp10(:,2),'g--',dm1ETFbp10(:,6),dm1ETFbp10(:,3),'g',dm1ETFbp10(:,6),dm1ETFbp10(:,2),'g--',dm2ETFbp10(:,6),dm2ETFbp10(:,3),'b',dm2ETFbp10(:,6),dm2ETFbp10(:,2),'b--',dm3ETFbp10(:,6),dm3ETFbp10(:,3),'r',dm3ETFbp10(:,6),dm3ETFbp10(:,2),'r--')

a=gca();
a.log_flags = "ln" ; 
mtlb_axis([0.01/100, 100/100, -6,5])

a=get("current_axes")
old_margins=a.margins
a.margins=old_margins

new_locations=a.x_ticks.locations

new_labels=string(new_locations*100)+"%"
a.x_ticks = tlist(["ticks", "locations", "labels"], new_locations, new_labels);


xs2eps(0,'FigQF1.eps')

//Figure 2: Performance (neg. excess returns=positive number :))
clf
plot(dm1ETFbp10(:,6),-dm1ETFbp10(:,5),'g--',dp2ETFbp10(:,6),-dp2ETFbp10(:,5),'g',dm2ETFbp10(:,6),-dm2ETFbp10(:,5),'b--',dp3ETFbp10(:,6),-dp3ETFbp10(:,5),'b',dm3ETFbp10(:,6),-dm3ETFbp10(:,5),'r--',dp4ETFbp10(:,6),-dp4ETFbp10(:,5),'r')
legend(["-1x";"+2x";"-2x";"+3x";"-3x";" +4x"],3)
//mtlb_axis([-0.03, 0, 0.005, 0.20]);
a=gca();
a.log_flags = "ll" ; 
mtlb_axis([0.1/100, 100/10^2,10^(-5)/100, 10/100])

a=get("current_axes")
old_margins=a.margins*0.8
a.margins=old_margins
//a.axes_reverse=["off","on","on"]

new_locations=a.x_ticks.locations
new_labels=string(new_locations*100)+"%"

a.x_ticks = tlist(["ticks", "locations", "labels"], new_locations, new_labels);


new_locations=a.y_ticks.locations
new_labels=string(new_locations*100)+"%"
a.y_ticks = tlist(["ticks", "locations", "labels"], new_locations,new_labels)

xs2eps(0,'FigQF2.eps')

    
    //Figure 3: gamma only 0.5 and 1, but large range of L. Vanishing Risk premium (kappa=0, therefore mu=0)
    //Computations are displayed on console gamma by gamma, but take a while.
    clf
    
    //Compute R^2 from trading boundaries for trading costs 10BP (could be changed, or added as input variable)
    function f=R2(gam,L)
    
      
        data=perf_trading(gam,0,L,0.001,0)
        pim=data(2)
        pip=data(3)
        
    f=pim*pip*(log(pip/pim)/(pip-pim))^2
    
    endfunction
    
    R2_vec=zeros(20,3)
    
    //Specify parameter ranges
    
    //multipliers
    L_vec=[-10:1:-5, -4.9:0.1:-0.1, 0.1:0.1:0.9, 1.1:0.1:4.9, 5:1:10];
    
    gam_1=0.5
    gam_2=1
    
    kappa=0
    bp=10
    
    for i=1:length(L_vec)
        
    L=L_vec(i)
    R2_vec(i,:)=[L, R2(gam_1,L),R2(gam_2,L)]//0.0506
    end
    
    
    i=i+1
    
    R2_vec(i,:)=[0,0,0]
    R2_vec(i+1,:)=[1,1,1]
    
    
    v=gsort(R2_vec,'lr','i')
    
    a=get("current_axes")
    old_margins=a.margins
    
    
    new_locations=[97.5,98,98.5,99.0,99.5,100]/100
    new_labels=string(new_locations*100)+"%"
    a.y_ticks = tlist(["ticks", "locations", "labels"], new_locations, new_labels);
    
    new_locations=[-10:1:10]
    new_labels=string(new_locations)
    a.x_ticks = tlist(["ticks", "locations", "labels"], new_locations, new_labels);
    
    
    plot(v(:,1),v(:,2),v(:,1),v(:,3))
    mtlb_axis([-10,10,0.975,1])
    
    
    xs2eps(0,'FigQF3.eps')

kappa1=0.04/0.16^2
kappa2=0.08/0.16^2



//Figures 4 and 5 
//Produce trading boundaries for riskpremia (uncomment to reproduce files with ending *1.5625 and ***3.125)
//gam=[900000:-100000:100000,99000:-1000:1000,900:-100:100,99.9:-0.1:0.5, 0.499:-0.01:0.09, 0.0899:-0.0001:0.0012,0.001201:-0.00001:0.000051];
//perf_trading(gam,0.08/0.16^2,3,0.001,1)

//perf_trading(gam,0.04/0.16^2,3,0.001,1)

//perf_trading(gam,0.08/0.16^2,-2,0.001,1)

//perf_trading(gam,0.04/0.16^2,-2,0.001,1)

clf
/////////10 BP

x1=fscanfMat('Bull30bp10kappa0');

x2=fscanfMat('Bull30bp10kappa1.5625');

x3=fscanfMat('Bull30bp10kappa3.125');


L=3
eps=1/1000
gam=x1(:,1)

N=length(gam)

//N=sum(gam<10)

gam=gam(1:N,:)

pim1=L-(3 ./gam/4 *L^2 *(L-1)^2 ).^(1/3)*0.001^(1/3)
pip1=L+(3 ./gam/4 *L^2 *(L-1)^2 ).^(1/3)*0.001^(1/3)

//subplot(2,1,1)
plot(x1(1:N,6),x1(1:N,2),'k')
plot(x1(1:N,6),x1(1:N,3),'k')
plot(x1(1:N,6),pim1,'r',x1(1:N,6),pip1,'r')



//plot(ziss,pim2,'g',ziss,pip2,'g')


gam=x2(:,1)
N=length(gam)
//N=sum(gam<10)
gam=gam(1:N,:)
plot(x2(1:N,6),x2,2),'b')
plot(x2(1:N,6),x2(1:N,3),'b')


gam=x3(:,1)
N=length(gam)
N=sum(gam<10)
gam=gam(1:N,:)
plot(x3(1:N,6),x3(1:N,2)*0.997,'g') 
plot(x3(1:N,6),x3(1:N,3),'g')

a=gca();
a.log_flags = "ln" ; 
mtlb_axis([0.001, 0.05, 2, 4])


a=get("current_axes")
new_locations=a.x_ticks.locations
new_labels=string(new_locations*100)+"%"
a.x_ticks = tlist(["ticks", "locations", "labels"], new_locations, new_labels);



xs2eps(0,'FigQF4.eps')



//For Figure 5, we need to also load the inverse funds
clf



y1=fscanfMat('Bear20bp10kappa0');
y2=fscanfMat('Bear20bp10kappa1.5625');
y3=fscanfMat('Bear20bp10kappa3.125');




epsi=10/10000
sigma=0.16
L=3
x1(:,7)=-3^(1/2)/12*sigma^3*L^2*(1-L)^2*epsi ./x1(:,6);

x2(:,7)=-3^(1/2)/12*sigma^3*L^2*(1-L)^2*epsi ./x2(:,6);

x3(:,7)=-3^(1/2)/12*sigma^3*L^2*(1-L)^2*epsi ./x3(:,6);

epsi=10/10000

L=3
y1(:,7)=-3^(1/2)/12*sigma^3*L^2*(1-L)^2*epsi ./y1(:,6);

y2(:,7)=-3^(1/2)/12*sigma^3*L^2*(1-L)^2*epsi ./y2(:,6);

y3(:,7)=-3^(1/2)/12*sigma^3*L^2*(1-L)^2*epsi ./y3(:,6);
//remove zero entries
y3=y3(4:length(y3(:,1)),:)

clf
plot(y1(:,6),-y1(:,5),'g-.',y2(:,6),-y2(:,5),'b-.',y3(:,6),-y3(:,5),'r-.',x1(:,6),-x1(:,5),'g',x2(:,6),-x2(:,5),'b',x3(:,6),-x3(:,5),'r')
plot(x1(:,6),-x1(:,7),'k')
a=gca();
a.log_flags = "ll" ; 
mtlb_axis([0.005, 10/100, 0.01/100, 0.01])

a=get("current_axes")
old_margins=a.margins
a.margins=old_margins*0.6

new_locations=[1 2 3 4 5 10]'/100
new_labels=string(new_locations*100)+"%"
a.x_ticks = tlist(["ticks", "locations", "labels"], new_locations, new_labels);

new_locations=[0.01 0.1 1]'/100
new_labels=string(new_locations*100)+"%"
a.y_ticks = tlist(["ticks", "locations", "labels"], new_locations, new_labels);

xs2eps(0,'FigQF5.eps')



