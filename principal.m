h = [0.3, 0.2, 0.1, 0.05, 0.03, 0.02];
nom_maillage1 = str2mat('geomCarre0.3.msh','geomCarre0.2.msh','geomCarre0.1.msh','geomCarre.05.msh','geomCarre.03.msh','geomCarre.02.msh');
h = [0.2, 0.12, 0.08, 0.03, 0.02, 0.015];
nom_maillage2 = str2mat('geomCarre_per0.2.msh','geomCarre_per.12.msh','geomCarre_per.08.msh','geomCarre_per.03.msh','geomCarre_per.02.msh','geomCarre_per.15.msh');
Neumann = 'non';
Dirichlet = 'oui';
Periodique = 'non';
erreur = zeros(6,2);
%[erreurL2,erreurL22,erreurH1,erreurH11] = principal_neumann('geomCarre.05.msh');
if strcmp(Neumann,'oui')
  for i = 1:length(h)
            nom=nom_maillage1(i,:);
            [erreurL2,erreurH1] = principal_neumann(nom);
            erreur(i,1)= erreurH1;
            erreur(i,2)= erreurL2;
  end
end
if strcmp(Periodique,'oui')
  for i = 1:length(h)
            nom=nom_maillage2(i,:);
            [erreurL2,erreurH1] = principal_periodique(nom);
            erreur(i,1)= erreurH1;
            erreur(i,2)= erreurL2;
  end
end
if strcmp(Dirichlet,'oui')
  for i = 1:length(h)
            nom=nom_maillage2(i,:);
            [erreurL2,erreurH1] = principal_dirichlet(nom);
            erreur(i,1)= erreurH1;
            erreur(i,2)= erreurL2;
  end
end
figure(1)
plot(-log(h),log(erreur(:,2)),'+-',-log(h),log(erreur(:,1)/(sqrt(2)*pi)),'o-');
xlabel('-logh'); ylabel('loge');
legend({'erreurL2', 'erreurH1'},'Location','southeast');