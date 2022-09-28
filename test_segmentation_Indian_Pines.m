
clear; close all; clc; 

load Indian_pines.mat

masque = double(masque);
I = double(I);

[Nx, Ny, Nbands] = size(I);
I = reshape(I, [Nx*Ny, Nbands]);

%---------------------------------
% R?duction de dimension par ACP
%---------------------------------
[coeff,score,latent] = acp(centree_reduire(I));
Nbands = 10; % Nombre de composantes principales conserv?es
I = score(:,1:Nbands);
I = reshape(I, [Nx, Ny, Nbands]);

% extraction des classes.
classes = unique(masque(masque~=0));
NbClasses = length(classes);

%-------------------------------------------------------------------------
% Estimation des param?tre des lois Gaussiennes multivari?es pour chacune
% des classes. 
%-------------------------------------------------------------------------
Mu_cl = cell(1, NbClasses); % Cell contenant les vecteurs moyens dans chaque classe.
Cov_cl = cell(1, NbClasses); % Cell contenant les matrices de covariance dans chaque classe.
for cl = 1:NbClasses
    [posx, posy] = find(masque==cl);
    I_cl = zeros(numel(posx),Nbands);
    for j = 1:numel(posx)
    I_cl(j,:) = squeeze(I(posx(j),posy(j),:));
    end
    % Calcul des moyennes dans chaque classe
    Mu_cl{cl} = mean(I_cl);
    % Calcul des matrices de covariance dans chaque classe
    Cov_cl{cl} = cov(I_cl);
end

%----------------------------------------
% Calcul des cartes de log-vraisemblance
%----------------------------------------
ln_p = zeros(size(I,1), size(I,2), NbClasses);
for cl=1:NbClasses
    % fonction calc_log_pdf_Gaussian à créer !
    ln_p(:,:,cl) = calc_log_pdf_Gaussian(I, Mu_cl{cl}, Cov_cl{cl});
end

%-----------------------------------------------
% Classification par maximum de vraisemblance
%-----------------------------------------------
disp('Classification par max de vraisemblance.')
[C,Id] = max(ln_p,[],3);

figure, imagesc(Id .* double(masque~=0)), caxis([0 NbClasses])
figure, imagesc(masque .* double(masque~=0)), caxis([0 NbClasses])

%-----------------------------
% R?gularisation Markovienne
%-----------------------------
disp('Segmentation Markovienne.')
EtiqInit = Id;
Beta = 1; % param?tre de r?gularisation
% fonction icm ? modifier pour prendre en entr?e les cartes de log-vraisemblance !
z = icm(ln_p, EtiqInit, NbClasses, Beta);
figure, imagesc(z .* double(masque~=0)), caxis([0 NbClasses])
