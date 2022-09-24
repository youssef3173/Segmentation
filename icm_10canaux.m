function EtiqResult = icm_10canaux(Im, EtiqInit, NbClasses, Mu, Sigma, Beta)
% Exécute une segmentation par champ de markov de l'image Im,
% en NbClasses classes dont les moyennes et les écart-types de niveaux de gris
% sont donnés par les tableaux Mu et Sigma
% EtiqInit, passée en paramètre, est l'image des étiquettes d'initialisation
% L'image des étiquettes en cours de modifications est Etiq
% L'image des étiquettes à renvoyer est EtiqResult

% L'algorithme utilisé est l'algorithme ICM (Iterated Conditional Modes),
% avec hypothèse de normalité (chaque classe suit une distribution  normale)

% L'énergie locale associée au site s sera notée Es. Elle sera composée de :
%		Es_rap : énergie de rappel aux données
%		Es_reg : énergie de régularisation
% La première porte sur la connaissance a priori des distributions des classes supposées gaussiennes
% La seconde permet de régulariser les régions. Elle s'appuiera sur un potentiel autologistique (de valeur +/- Beta)


Im=double(Im);
[Nb_lignes,Nb_colonnes]=size(Im);
Etiq = EtiqInit;

% Début du processus itératif
k=0;
fin=0;

while fin==0
    % Nombre de modifications sur un balayage complet
    nb_modifs=0;
    % Parcours de l'image (les bords ne seront pas pris en compte dans un premier temps)
    for i=2:Nb_lignes-1
        for j=2:Nb_colonnes-1
            
            % Calcul de l'énergie locale Es pour toutes les étiquettes lambda possibles
            
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
            
            % Recherche de la meilleure étiquette

            
            min = Etot(1);
            idx = 1;
            for c = 2:NbClasses
                if( Etot(c) < min )
                   idx = c;
                   min = Etot(c);
                end
            end
            
            % Affectation de la meilleure étiquette

            if (  Etiq( i, j) ~= idx )
                nb_modifs = nb_modifs + 1;
            end
            Etiq( i, j) = idx;
            
        end
    end
    
    k=k+1;
    if (nb_modifs==0)
        fin=1; % On met fin à la boucle s'il n'y a plus de changement.
        
        figure,
        subplot 121, imshow( EtiqInit , colormap( jet(4))), title('Initial');
        subplot 122, imshow( Etiq, colormap( jet(4))), title('Resultat');
    end
    fprintf('Iteration : %6d      -      Nb de modifications : %6d\n',k,nb_modifs);
    
    % Affichage de l'image intermédiaire des étiquettes
    % A compléter...
end
EtiqResult=Etiq;

end


