% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Neumann sur le maillage nom_maillage.msh
%
% | -Delta  u + u= f,   dans \Omega
% |         du/dn = 0,   sur le bord
%
% =====================================================

function [errL2,errH1]=principal_neumann(nom_maillage)
% lecture du maillage et affichage
% ---------------------------------
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
  
   Kel=matK_elem1(S1, S2, S3);
           
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

% inversion
% ----------
UU = (MM+KK)\LL;

% visualisation
% -------------
validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
%UU_exact = cos(pi*Coorneu(:,1)).*cos(2*pi*Coorneu(:,2));
UU_exact = cos(2*pi*Coorneu(:,1));
% Calcul de l erreur L2
% A COMPLETER
errL2 = sqrt((UU_exact-UU)' * MM * (UU_exact-UU));
% Calcul de l erreur H1
errH1 = sqrt((UU_exact-UU)' * KK * (UU_exact-UU));
% attention de bien changer le terme source (dans FF)
end
affiche(UU_exact, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));
affiche(UU, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
