# Projet algèbre linéaire distribuée
## Sujet 3 : Deeper look on Hypre BoomerAMG

### Pour compiler le projet :
- Cloner le [repository Hypre](https://github.com/hypre-space/hypre)
- Compiler hypre en suivant [la notice d'installation](https://github.com/hypre-space/hypre/blob/master/INSTALL.md) 
    - Penser à export CC, CXX et FC avec les compilateurs mpi (comme dans [cet exemple d'erreur](https://stackoverflow.com/questions/26920083/fatal-error-mpi-h-no-such-file-or-directory-include-mpi-h))
- Positionner la variable `HYPRE_DIR` dans le fichier `Makefile`
- Si vous êtes sur le cluster du CEMEF, il faudra peut être changer les modules correspondants avec le script `load_modules`
- Lancer la compilation avec `make tp`


### Paramètre `-file`
A l'éxécution, vous pouvez utiliser l'option `-file chemin/vers/matrixmarket.mtx` pour utiliser le solveur sur une matrice réelle. 
Vous trouverez quelques matrices dans le répertoire `matrices/`.

Sinon, n'hésitez pas à télécharger une matrice symmétrique définie positive au format CSR sur [ce site](https://sparse.tamu.edu/).
