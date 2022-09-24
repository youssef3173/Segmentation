function EtiqResult = icm_10canaux(Im, EtiqInit, NbClasses, Mu, Sigma, Beta)
% Ex�cute une segmentation par champ de markov de l'image Im,
% en NbClasses classes dont les moyennes et les �cart-types de niveaux de gris
% sont donn�s par les tableaux Mu et Sigma
% EtiqInit, pass�e en param�tre, est l'image des �tiquettes d'initialisation
% L'image des �tiquettes en cours de modifications est Etiq
% L'image des �tiquettes � renvoyer est EtiqResult

% L'algorithme utilis� est l'algorithme ICM (Iterated Conditional Modes),
% avec hypoth�se de normalit� (chaque classe suit une distribution  normale)

% L'�nergie locale associ�e au site s sera not�e Es. Elle sera compos�e de :
%		Es_rap : �nergie de rappel aux donn�es
%		Es_reg : �nergie de r�gularisation
% La premi�re porte sur la connaissance a priori des distributions des classes suppos�es gaussiennes
% La seconde permet de r�gulariser les r�gions. Elle s'appuiera sur un potentiel autologistique (de valeur +/- Beta)


Im=double(Im);
[Nb_lignes,Nb_colonnes]=size(Im);
Etiq = EtiqInit;

% D�but du processus it�ratif
k=0;
fin=0;

while fin==0
    % Nombre de modifications sur un balayage complet
    nb_modifs=0;
    % Parcours de l'image (les bords ne seront pas pris en compte dans un premier temps)
    for i=2:Nb_lignes-1
        for j=2:Nb_colonnes-1
            
            % Calcul de l'�nergie locale Es pour toutes les �tiquettes lambda possibles
            
            E1 = zeros( 1, NbClasses);
            Etot = zeros( 1, NbClasses);
            for c = 1:NbClasses
                U1 = 0; 
                l1 = c;  % ettiquite de test
                
                % test sur les 4 voisins :
                l2 = [ Etiq( i + 1, j), Etiq( i - 1, j), Etiq( i , j + 1), Etiq( i , j - 1) ];
                % test sur les 8 voisins :
                % l2 = [ Etiq( i + 1, j), Etiq( i - 1, j), Etiq( i , j + 1), Etiq( i , j - 1), ...
                %        Etiq( i + 1, j + 1), Etiq( i - 1, j + 1), Etiq( i + 1, j - 1), Etiq( i - 1, j - 1) ];
                
                for elt = 1:length(l2)
                    if ( l1 == l2(elt))
                        U1 = U1 - Beta; 
                    else 
                        U1 = U1 + Beta;
                    end  
                end
                
                E1(c) = U1;
                
                vect = zeros( 1, 10);
                for count = 1:10
                   vect( count) = Im(x,y,i); 
                end
                E2 = log( (2*pi)^(d/2)*det(sigma) ) + 0.5*(vect - Mu)*pinv( Sigma)*(vect - Mu)';

                % Energie totale :
                Etot(c) = E2 + E1(c); 
                
            end
            
            % Recherche de la meilleure �tiquette

            
            min = Etot(1);
            idx = 1;
            for c = 2:NbClasses
                if( Etot(c) < min )
                   idx = c;
                   min = Etot(c);
                end
            end
            
            % Affectation de la meilleure �tiquette

            if (  Etiq( i, j) ~= idx )
                nb_modifs = nb_modifs + 1;
            end
            Etiq( i, j) = idx;
            
        end
    end
    
    k=k+1;
    if (nb_modifs==0)
        fin=1; % On met fin � la boucle s'il n'y a plus de changement.
        
        figure,
        subplot 121, imshow( EtiqInit , colormap( jet(4))), title('Initial');
        subplot 122, imshow( Etiq, colormap( jet(4))), title('Resultat');
    end
    fprintf('Iteration : %6d      -      Nb de modifications : %6d\n',k,nb_modifs);
    
    % Affichage de l'image interm�diaire des �tiquettes
    % A compl�ter...
end
EtiqResult=Etiq;

end


