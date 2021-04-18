% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions periodiques
% sur le maillage nom_maillage.msh
%
% | -div(A grad u ) u + u= f,   dans \Omega
% |         u periodique,   sur le bord
%
% =====================================================


% lecture du maillage et affichage
% ---------------------------------
function [errL2,errH1]=principal_periodique(nom_maillage)
%nom_maillage = 'geomCarre_per.05.msh';
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
    
  % On fait l'assemblage de la matrice globale et du second membre
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

% Projection sur l espace V_p
% ———————————————————
% matrice de projection 
[K]=find(Refneu == 0);
[I]=find(Refneu == 5);
[J]=find(Refneu == 1);
[JJ]=find(Refneu == 3);
[Q]=find(Refneu == 2);
[QQ]=find(Refneu == 4);
n = (Nbpt - length(K))/2;
PP = sparse(Nbpt - n - 1, Nbpt);
for i = 1:length(I)
      PP(1,I(i))= 1;
end
for  i = 1:length(J)
    PP(i+1,J(i))= 1;
    for j = 1:length(J)
        if abs(Coorneu(JJ(j),1) - Coorneu(J(i),1)) < 0.005
            PP(i+1,JJ(j))= 1;
        end
    end   
end

for  i = 1:length(Q)
    PP(i+1+length(J),Q(i))= 1;
    for j = 1:length(J)
        if abs(Coorneu(QQ(j),2) - Coorneu(Q(i),2))< 0.005
            PP(i+1+length(J),QQ(j))= 1;
        end
    end   
end

for i = n: Nbpt - n - 1
     PP(i,K(i-n+1)) = 1;
end
AA = MM+KK;
AAp = PP*AA*transpose(PP);
LLp = PP*LL;

% inversion
% ----------
UUp = AAp\LLp;


% Expression de la solution dans toute la base
% ———————
UU = transpose(PP)*UUp;

% visualisation
% -------------
UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));
affiche(UU, Numtri, Coorneu, sprintf('Periodique - %s', nom_maillage));
affiche(UU_exact, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));
validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
% Calcul de l erreur L2
errL2 = sqrt((UU_exact-UU)' * MM * (UU_exact-UU));
% Calcul de l erreur H1
errH1 = sqrt((UU_exact-UU)' * KK * (UU_exact-UU));
%bien changer le terme source (dans FF)
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

