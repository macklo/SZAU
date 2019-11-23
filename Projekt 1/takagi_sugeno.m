function [a, c, h2r0] = takagi_sugeno(il, draw)
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
    V1r0 = A1*(alfa2/alfa1)^2*h2r0;
    V2r0 = C2.*h2r0.*h2r0;
    F1r0 = alfa1.*h1r0.^0.5-FD;
%     simgmf(ymin:0.1:ymax , [-3 c(1)])
%    dsigmf(ymin:0.1:ymax, [3 c(r-1) -3 c(r)]
    
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
        ylabel('przynale¿noœæ')
        title('Funkcje przynale¿noœci regulatorów lokalnych w regulacji rozmytej wzglêdem wartoœci wyjœcia')
    end
    F10 = 90;
    FD0 = 30;
    h20 = 36;
    h10 = (alfa2/alfa1)^2*h20;
    V10 = A1*h10;
    V20 = C2*h20*h20;

    if draw
        figure
        title('Przebiegi wyjœcia dla skoku wartoœci sterowania')
        xlabel('czas[t]')
        ylabel('h2[cm]')
        hold on
    end
    for i = 180:-18:0
        
        F1in = F10 * ones(1,n);
        F1in(start:n) = i;
        F1 = F10 * ones(1,n);
        FD = FD0 * ones(1,n);
        h1 = h10 * ones(2+il+1,n);
        h2 = h20 * ones(2+il+1,n);
        V1 = A1*h10 * ones(2+il+1,n);
        V2 = C2*h20*h20 * ones(2+il+1,n);
        w = ones(1,il);
        for t = tau+1 : n
            F1(t) = F1in(t-tau);
            
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
                
%                 V1(2+r,t) = V1(2+il+1,t-1) + (F1(t-1) - Fr0(r)) + (FD(t-1) - FD0) - a1/2*hr0(r)^-0.5 * (h1(2+il+1,t-1) - hr0(r));
%                 V2(2+r,t) = V2(2+il+1,t-1) + a1/2*hr0(r)^-0.5 * (h1(2+il+1,t-1) - hr0(r)) - a1/2*hr0(r)^-0.5 * (h2(2+il+1,t-1) - hr0(r));
%                 h1(2+r,t) = hr0(r) + 1/2*(C1*Vr0(r))^-0.5 * (V1(2+r,t) - Vr0(r));
%                 h2(2+r,t) = hr0(r) + 1/2*(C2*Vr0(r))^-0.5 * (V2(2+r,t) - Vr0(r));
                
                
                %obliczanie wag
               
                if r == 1
                    w(r) = sigmf(h2(2+il+1,t-1) , [-a c(1)]);

                    
                elseif r == il
                    w(r) = sigmf(h2(2+il+1,t-1) , [a c(il-1)]);
                else
                    w(r) = dsigmf(h2(2+il+1,t-1), [a c(r-1) a c(r)]);
                end
            end


            
            %obliczanie wyjœæ modelu rozmytego
            h2(2+il+1,t) = w*h2(3:2+il, t)/sum(w);
            h1(2+il+1,t) = w*h1(3:2+il, t)/sum(w);
            V1(2+il+1,t) = w*V1(3:2+il, t)/sum(w);
            V2(2+il+1,t) = w*V2(3:2+il, t)/sum(w);
        end
        if draw
            plot(start:n,h2(1,start:end),'g')
            plot(start:n,h2(2,start:end),'r')
            plot(start:n,h2(2+il+1,start:end),'b')
        end
    end
end
