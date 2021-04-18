% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Dirichlet sur le maillage nom_maillage.msh
%
% | -\Delta u + u= f,   dans \Omega
% |         u = 0,   sur le bord
%
% =====================================================


% lecture du maillage et affichage
% ---------------------------------
function [errL2,errH1]=principal_dirichlet(nom_maillage)
%nom_maillage = 'geomCarre0.1.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  % A COMPLETER
  S1=Coorneu(Numtri(l,1),:);
  S2=Coorneu(Numtri(l,2),:);
  S3=Coorneu(Numtri(l,3),:);
  % calcul des matrices elementaires du triangle l 
  
   Kel=matK_elem(S1, S2, S3);
           
   Mel=matM_elem(S1, S2, S3);
    
  % On fait l'assemmblage de la matrice globale et du second membre
  % A COMPLETER
  for i = 1:3
        for j = 1:3
            MM(Numtri(l,i),Numtri(l,j)) = MM(Numtri(l,i),Numtri(l,j)) + Mel(i,j);
            KK(Numtri(l,i),Numtri(l,j)) = KK(Numtri(l,i),Numtri(l,j)) + Kel(i,j);
        end % for j
  end  % for i

end % for l

% Calcul du second membre L
% -------------------------
	% A COMPLETER
	% utiliser la routine f.m
FF = f(Coorneu(:,1),Coorneu(:,2));
LL = MM*FF;

% Projection sur l espace V_0
% ———————————————————
% matrice de projection 

[J]=find(Refneu == 0);
PP = sparse(length(J),Nbpt);
for i = 1:length(J)
     PP(i,J(i)) = 1;
end
AA = MM+KK;
AA0 = PP*AA*transpose(PP);
LL0 = PP*LL;

% inversion
% ----------
UU0 = AA0\LL0;

% Expression de la solution dans toute la base
% ———————
UU = transpose(PP)*UU0;

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));

validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));
% Calcul de l erreur L2
errL2 = sqrt((UU_exact-UU)' * MM * (UU_exact-UU));
% Calcul de l erreur H1
errH1 = sqrt((UU_exact-UU)' * KK * (UU_exact-UU));
% attention de bien changer le terme source (dans FF)
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

