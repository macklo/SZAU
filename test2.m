








addpath("./classes")
addpath("./abstraction")
il =2;
draw =1;
n = 4000;
A1 = 540;
C2 = 0.85;
alfa1 = 26;
alfa2 = 20;
tau = 100;
start = 1;
FD = 30;

ymax = 110.2;
ymin = 2.25;
dy = (ymax-ymin)/il;
a = 0.5;
c = ymin+dy:dy:ymax-dy;
h2r0 = ones(1,il);
h2r0(1) = (c(1)+ymin)/2-1;
h2r0(il) = min((ymax+c(il-1))/2+1, ymax);
h1r0 = (alfa2/alfa1)^2*h2r0;
if il > 2
    h2r0(2:il-1) = (c(2:il-1)+c(1:il-2))./2;
    h1r0(2:il-1) = (alfa2/alfa1)^2*h2r0(2:il-1);
    
end

for i= 1:il
    workpoints(i) = calculate_workpoint(h2r0(i));
%     lintanks(i) = LinearTankSystem(workpoints(1));
%     lintanks(i).resetToWorkPoint(workpoints(1));
    
    
end

lintanks1 = LinearTankSystem(workpoints(1))
lintanks2 = LinearTankSystem(workpoints(2))
lintanks1.resetToWorkPoint(workpoints(1));
lintanks2.resetToWorkPoint(workpoints(2));

V1r0 = A1*(alfa2/alfa1)^2*h2r0;
V2r0 = C2.*h2r0.*h2r0;
F1r0 = alfa1.*h1r0.^0.5-FD;

if draw
    figure
    hold on
    for r = 1:il
        if r == 1
            plot(ymin:0.1:ymax, sigmf(ymin:0.1:ymax , [-a c(1)]))
        elseif r == il
            plot(ymin:0.1:ymax,sigmf(ymin:0.1:ymax , [a c(il-1)]))
        else
            plot(ymin:0.1:ymax,dsigmf(ymin:0.1:ymax, [a c(r-1) a c(r)]))
        end
    end
    plot(h2r0, ones(1,il), 'ko')
    xlabel('h2')
    ylabel('przynale?no??')
    title('Funkcje przynale?no?ci regulator?w lokalnych w regulacji rozmytej wzgl?dem warto?ci wyj?cia')
end
F10 = 90;
FD0 = 30;
h20 = 36;
h10 = (alfa2/alfa1)^2*h20;
V10 = A1*h10;
V20 = C2*h20*h20;

if draw
    figure
    title('Przebiegi wyj?cia dla skoku warto?ci sterowania')
    xlabel('czas[t]')
    ylabel('h2[cm]')
    hold on
end

workpoint = calculate_workpoint(36);
tank       =TankSystem(workpoint);
tank.resetToWorkPoint(workpoint);

lintank  = LinearTankSystem(workpoint);
lintank.resetToWorkPoint(workpoint);
u11 = workpoint.u*ones(1,n);


for i = 135:135
    
    u11(start:end)=i;
    F1in = F10 * ones(1,n);
    F1in(start:n) = i;
    F1 = F10 * ones(1,n);
    F1(101:end) =i ;
    FD = FD0 * ones(1,n);
    h1 = h10 * ones(2+il+1,n);
    h2 = h20 * ones(2+il+1,n);
    V1 = A1*h10 * ones(2+il+1,n);
    V2 = C2*h20*h20 * ones(2+il+1,n);
    w = ones(1,il);
    dynout = workpoint.y*ones(1,n);
    linout = workpoint.y*ones(1,n);
    fuzz(1,:) = workpoint.y*ones(1,n);
    fuzz(2,:) = workpoint.y*ones(1,n);
    fuzzout(1,:) = workpoint.y*ones(1,n);
    for t = 2 : n
        
        %model dynamiczny
        tank.setControl(u11(t));
        dynout(t)   = tank.getOutput();
        tank.nextIteration();
        
        %model liniowy
        lintank.setControl(u11(t));
        linout(t)   = lintank.getOutput();
        lintank.nextIteration();
        
        %modele liniowe modelu rozmytego
        lintanks1.setControl(u11(t));
        lintanks2.setControl(u11(t));
        fuzz(1,t)   = lintanks1.getOutput();
        fuzz(2,t)   = lintanks2.getOutput();
        lintanks1.nextIteration();
        lintanks2.nextIteration();
        
        
        
        wk(1) = sigmf(fuzz(1,t) , [-a c(1)]);
        wk(2) = sigmf(fuzz(2,t) , [a c(il-1)]);
        fuzzout(t) =  wk*fuzz(:, t)/sum(wk);
        
        
        %model dynamiczny
        V1(1,t) = V1(1,t-1) + F1(t-1)+ FD(t-1) - alfa1*h1(1,t-1)^0.5;
        V2(1,t) = V2(1,t-1) + alfa1*h1(1,t-1)^0.5 - alfa2*h2(1,t-1)^0.5;
        h1(1,t) = V1(1,t)/A1;
        h2(1,t) = (V2(1,t)/C2)^0.5;
        
        %model liniowy
        V1(2,t) = V1(2,t-1) + (F1(t-1)-F10) + (FD(t-1)-FD0)   - alfa1/2*h10^-0.5 * (h1(2,t-1) - h10);
        V2(2,t) = V2(2,t-1)   + alfa1/2*h10^-0.5 * (h1(2,t-1) - h10)    - alfa2/2*h20^-0.5 * (h2(2,t-1) - h20);
        h1(2,t) = (V1(2,t))/A1;
        h2(2,t) = h20 + 1/2*(C2*V20)^-0.5 * (V2(2,t) - V20);
        
        for r = 1:il
            %modele liniowe modelu rozmytego
            V1(2+r,t) = V1(2+il+1,t-1) + (F1(t-1)-F1r0(r)) + (FD(t-1)-FD0)   - alfa1/2*h1r0(r)^-0.5 * (h1(2+il+1,t-1) - h1r0(r));
            V2(2+r,t) = V2(2+il+1,t-1) + alfa1/2*h1r0(r)^-0.5 * (h1(2+il+1,t-1) - h1r0(r))    - alfa2/2*h2r0(r)^-0.5 * (h2(2+il+1,t-1) - h2r0(r));
            h1(2+r,t) = h1r0(r) + (V1(2+r,t)- V1r0(r))/A1;
            h2(2+r,t) = h2r0(r) + 1/2*(C2*V2r0(r))^-0.5 * (V2(2+r,t) - V2r0(r));
            
            
            
            %obliczanie wag
            
            if r == 1
                w(r) = sigmf(h2(2+il+1,t-1) , [-a c(1)]);
                
                
            elseif r == il
                w(r) = sigmf(h2(2+il+1,t-1) , [a c(il-1)]);
            else
                w(r) = dsigmf(h2(2+il+1,t-1), [a c(r-1) a c(r)]);
            end
        end
        
        
        
        %obliczanie wyj?? modelu rozmytego
        h2(2+il+1,t) = w*h2(3:2+il, t)/sum(w);
        h1(2+il+1,t) = w*h1(3:2+il, t)/sum(w);
        V1(2+il+1,t) = w*V1(3:2+il, t)/sum(w);
        V2(2+il+1,t) = w*V2(3:2+il, t)/sum(w);
    end
    
end
if draw
    plot(h2(1,start:end),'g')
    plot(h2(2,start:end),'r')
    plot(h2(2+il+1,start:end),'b')
    plot(dynout(start:end), 'y')
    plot(linout(start:end), 'm')
    plot(fuzzout(1,(start:end)), 'k')
%     plot(fuzz(1,(start:end)), 'k')
%     plot(fuzz(2,(start:end)), 'k')
    
end