%% Created by P. Heins 27.05.2015, University of Sheffield

%% This program takes in output-feedback LTI controller matrices from multiple
%% controllers and prints them to files in binary format, ready to be read 
%% by the controller class in Channelflow.
%%
%% Controllers should be of form:
%% \dot{x} = Ax + By
%%       u = Cx + Dy
%% x - controller state vector, y - measurements from simulation
%% u - control inputs to simulation
%%
%% Each controller should be named K(m/p)kx(m/p)kz, where m/p
%% denotes minus/positive. E.g. controller for kx=1,kz=1 will
%% be named Kp1p1, controller for kx=-2,kz=3 will be named
%% Km2p3 etc.



% Parameters
minkx = 0;
maxkx = 0;
minkz = 0;
maxkz = 10;
Ny_Sim = 151;   % Simulation wall-normal resolution
Re = 2230;

Lx = 2*pi;
Lz = 2*pi;


Mx = (maxkx-minkx)+1;
Mz = (maxkz-minkz)+1;

StateArray = zeros(Mx*Mz,3);

mkdir('Mult_Control_Mat');



count = 0;
for kx = minkx:maxkx
    for kz = minkz:maxkz
        
        count = count + 1;
        
        if (kx==0 && kz==0)           
            
            continue
        end
        
               
        
        
        disp(['Kx=' num2str(kx) 'Kz=' num2str(kz)]);
  
        
        % get variable names
        if (kx <0)           
            K = sprintf('Km%ip%i',abs(kx),kz);           
        else 
            K = sprintf('Kp%ip%i',kx,kz);    
        end
        
     % Assign variable   
     K = eval(K);
     
     
     alpha = 2*pi*kx/Lx;
     beta = 2*pi*kz/Lz;
          
     Cs = getC(Ny_Sim,Re,alpha,beta);  % Sensor matrix - please create own function



     Amat = K.a;    
     Bmat = K.b;
     Cmat = K.c;
     Dmat = K.d;
       
     % Generate filenames
     fida = fopen(['Mult_Control_Mat/Ka_mat_kx' num2str(kx) 'kz' num2str(kz) '.bin'],'w');
     fidb = fopen(['Mult_Control_Mat/Kb_mat_kx' num2str(kx) 'kz' num2str(kz) '.bin'],'w');
     fidc = fopen(['Mult_Control_Mat/Kc_mat_kx' num2str(kx) 'kz' num2str(kz) '.bin'],'w');
     fidd = fopen(['Mult_Control_Mat/Kd_mat_kx' num2str(kx) 'kz' num2str(kz) '.bin'],'w');
     fidcs = fopen(['Mult_Control_Mat/Cs_mat_kx' num2str(kx) 'kz' num2str(kz) '.bin'],'w');

     % Write matrix dimensions at top of file in binary i.e. [row;col]   
     fwrite(fida,[size(Amat,1) size(Amat,2)],'double'); 
     fwrite(fidb,[size(Bmat,1) size(Bmat,2)],'double'); 
     fwrite(fidc,[size(Cmat,1) size(Cmat,2)],'double'); 
     fwrite(fidd,[size(Dmat,1) size(Dmat,2)],'double'); 
     fwrite(fidcs,[size(Cs,1) size(Cs,2)],'double'); 

     % Print matrices to file in binary
     for rw = 1:size(Amat,1)
         for col = 1:size(Amat,2)
             fwrite(fida,Amat(rw,col),'double');             
         end
     end
     
     for rw = 1:size(Bmat,1)
         for col = 1:size(Bmat,2)
             fwrite(fidb,Bmat(rw,col),'double');             
         end
     end
     
     for rw = 1:size(Cmat,1)
         for col = 1:size(Cmat,2)
             fwrite(fidc,Cmat(rw,col),'double');             
         end
     end
     
     for rw = 1:size(Dmat,1)
         for col = 1:size(Dmat,2)
             fwrite(fidd,Dmat(rw,col),'double');             
         end
     end
     
     for rw = 1:size(Cs,1)
         for col = 1:size(Cs,2)
             fwrite(fidcs,Cs(rw,col),'double');             
         end
     end
        
     
     fclose(fida);
     fclose(fidb);
     fclose(fidc);
     fclose(fidd);
     fclose(fidcs);
     
     
     % Calculate appropriate mx,mz mode numbers for channelflow
     if (kx<0)
         mx = kx+Mx;
     else
         mx = kx;
     end
     
     mz = kz;
     
     % records number of states for each mx,mz pair
     NStates = size(Amat,1);
     
     StateArray(count,1) = mx;
     StateArray(count,2) = mz;
     StateArray(count,3) = NStates;

         
     
     
     
    end
end

StateArray = sortrows(StateArray);

fidS = fopen(['Mult_Control_Mat/StateInfo.bin'],'w');
fwrite(fida,[size(StateArray,1) size(StateArray,2)],'double'); 

for rw = 1:size(StateArray,1)
         for col = 1:size(StateArray,2)
             fwrite(fidS,StateArray(rw,col),'double');             
         end
end

fclose(fidS);

fidS2 = fopen(['Mult_Control_Mat/Controller_Info.asc'],'w');
fprintf(fidS2,'Re=%g\n',Re);
fprintf(fidS2,'Tau_phi=%g\n',tau_phi);
fprintf(fidS2,'Lx=%g\n',Lx);
fprintf(fidS2,'Lz=%g\n',Lz);
fprintf(fidS2,'Kx range:%i-%i\n',minkx,maxkx);
fprintf(fidS2,'Kz range:%i-%i\n',minkz,maxkz);
fprintf(fidS2,'Ny_Sim=%g\n',Ny_Sim);
fprintf(fidS2,'Ny_Controller=%g\n',Ny);
fprintf(fidS2,'Continuous/Discrete: Continuous\n');
fprintf(fidS2,'Size(Ka)=%ix%i\n',size(Amat,1),size(Amat,2));
fprintf(fidS2,'Size(Kb)=%ix%i\n',size(Bmat,1),size(Bmat,2));
fprintf(fidS2,'Size(Kc)=%ix%i\n',size(Cmat,1),size(Cmat,2));
fprintf(fidS2,'Size(Kd)=%ix%i\n',size(Dmat,1),size(Dmat,2));
fprintf(fidS2,'Size(Cs)=%ix%i\n',size(Cs,1),size(Cs,2));
fclose(fidS2);
        
        
        

