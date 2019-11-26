clear
clc

addpath("./classes")
addpath("./abstraction")


workpoint = calculate_workpoint(36)

sim_length = 3000;
jumpK      = 100;
fd         = 30;

jumps  = -1:0.2:1;
uJumps = workpoint.u + jumps * workpoint.u;
dJumps = fd + jumps*fd;

y       = cell(size(jumps));

legends = cell(size(jumps));

ystat    = zeros(size(jumps));


tanks       = TankSystem(workpoint);


for i = 1:size(uJumps, 2)
    uJump = uJumps(i)
    dJump = dJumps(i)
   
    legends{i} = "F_{D} = " + num2str(dJump);
    
    tanks.resetToWorkPoint(workpoint);
    
    
    u = workpoint.u.*ones(1, sim_length);
   
    u(1, jumpK:end) = 90;
    tanks.setDisturbance(dJump);
    
    y{i}    = workpoint.y.*ones(1, sim_length);
    
    
    for k = 1:1:sim_length
        y{i}(k)    = tanks.getOutput();
        
        tanks.setControl(u(k));
        
        tanks.nextIteration();
        
    end
end


figure(1)
grid on
hold on
for i = 1:size(uJumps, 2)
    stairs(y{i})
end

set(gca, 'ColorOrderIndex', 1)


title ("Przebiegi wyjœcia modeli dla ró¿nych zak³óceñ")
xlabel("t [s]")
ylabel("h_2 [cm]")
legend(legends, 'Location', 'EastOutside')
%     print(figure(1),'rys3', '-dmeta' , '-r500' )


for i = 1:size(uJumps, 2)
    uJump = uJumps(i)

    legends{i} = "F_{1in} = " + num2str(uJump);
  
    
    tanks.resetToWorkPoint(workpoint);
    
    
    u = workpoint.u.*ones(1, sim_length);
    u(1, jumpK:end) = uJump;
    
         tanks.setDisturbance(30);
    
    y{i}    = workpoint.y.*ones(1, sim_length);
    
    
    for k = 1:1:sim_length
        y{i}(k)    = tanks.getOutput();
        
        tanks.setControl(u(k));
        
        tanks.nextIteration();
        
    end
end

figure(2)
	grid on
	hold on
	for i = 1:size(uJumps, 2)
		stairs(y{i})
	end
	
	set(gca, 'ColorOrderIndex', 1)
	
 	title ("Przebiegi wyjœcia modeli dla ró¿nych skoków sterowania")
	xlabel("t [s]")
	ylabel("h_2 [cm]")
 	legend(legends, 'Location', 'EastOutside')
%     print(figure(2),'rys3', '-dmeta' , '-r500' )
