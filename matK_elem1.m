function [Kel] = matK_elem1(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul la matrices de raideur elementaire en P1 lagrange
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exacte (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par 
%          les coordonnees barycentriques 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);
Bl = [x2-x1 x3-x1; y2-y1 y3-y1]; 
Sl = [x1 ; y1];
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
BL = 1/abs(D)*[y3-y1 y1-y2; x1-x3 x2-x1];
% les 3 normales a l'arete opposees (de la longueur de 'arete)
norm = [-1 -1;1 0;0 1];
M = zeros(4,2);
M(1,:) = Bl*[1/3;1/3]+Sl;
M(2,:) = Bl*[1/5;1/5]+Sl;
M(3,:) = Bl*[1/5;3/5]+Sl;
M(4,:) = Bl*[3/5;1/5]+Sl;
W = [-9/32,25/96,25/96,25/96];
% D est, au signe pres, deux fois l'aire du triangle
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end
% calcul de la matrice de raideur
% -------------------------------
Kel = zeros(3,3);
for i=1:3
  for j=1:3
	% A COMPLETER
    for k =1:4
     Kel(i,j) = Kel(i,j)+W(k)*(A(M(k,1),M(k,2))*(BL*norm(i,:)'))'*(BL*norm(j,:)')*abs(D); 
    end
  end % j
end % i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
