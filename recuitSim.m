function EtiqResult=recuitSim(Im, EtiqInit, NbClasses, Mu, Sigma, beta)



Im=double(Im);
[Nb_lignes,Nb_colonnes]=size(Im);
Etiq=EtiqInit;

% Début du processus itératif
iter=0;
fin=0;


while fin==0
    % Nombre de modifications sur un balayage complet
    nb_modifs=0;
    % Parcours de l'image (les bords ne seront pas pris en compte dans un premier temps)
    for i=2:Nb_lignes-1
        for j=2:Nb_colonnes-1 
            T = 5; 
            cur = Etiq(i,j);
            rd = randi([1 NbClasses]); 
            while (cur~=rd)
                rd = randi([1 NbClasses]);
            end
            
            Ex = log(Sigma(rd)*sqrt(2*pi)) + (Im(i,j)-Mu(rd))^2/(2*Sigma(rd)^2);
            Es = log(Sigma(cur)*sqrt(2*pi)) + (Im(i,j)-Mu(cur))^2/(2*Sigma(cur)^2);
            
            if Etiq(i-1,j) == cur
                Es = Es - beta; 
            else 
                Es = Es + beta; 
            end

            if Etiq(i+1,j) == cur
                Es = Es - beta; 
            else 
                Es = Es + beta; 
            end

            if Etiq(i,j-1) == cur
                Es = Es - beta; 
            else 
                Es = Es + beta; 
            end

            if Etiq(i,j+1) == cur
                Es = Es - beta; 
            else 
                Es = Es + beta; 
            end
            
            
            if Etiq(i-1,j) == rd
                Ex = Ex - beta; 
            else 
                Ex = Ex + beta; 
            end

            if Etiq(i+1,j) == rd
                Ex = Ex - beta; 
            else 
                Ex = Ex + beta; 
            end

            if Etiq(i,j-1) == rd
                Ex = Ex - beta; 
            else 
                Ex = Ex + beta; 
            end

            if Etiq(i,j+1) == rd
                Ex = Ex - beta; 
            else 
                Ex = Ex + beta; 
            end
            
            DeltaE = abs(Es-Ex); 
            if (Es>Ex)
                Etiq(i,j)=rd; 
                nb_modifs = nb_modifs+1; 
            else 
                if (exp(-DeltaE/T)>rand)
                    Etiq(i,j)=rd;
                    nb_modifs = nb_modifs+1; 
                end
            end
 
        end
    end
    T = 0.8*T; 
end
    
    iter=iter+1;
    if (nb_modifs==0)
        fin=1; % On met fin à la boucle s'il n'y a plus de changement.
    end
    fprintf('Iteration : %6d      -      Nb de modifications : %6d\n',iter,nb_modifs);
    
    % Affichage de l'image intermédiaire des étiquettes
    imshow(Etiq, colormap(jet(4))); 
    EtiqResult=Etiq; 
end




