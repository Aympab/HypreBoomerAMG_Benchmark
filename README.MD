# Projet algèbre linéaire distribuée
## Sujet 3 : Deeper look on Hypre BoomerAMG

Pour compiler le projet :
- Cloner le [repository Hypre](https://github.com/hypre-space/hypre)
- Compiler hypre en suivant [la notice d'installation](https://github.com/hypre-space/hypre/blob/master/INSTALL.md) 
    - Penser à export CC, CXX et FC avec les compilateurs mpi (comme dans [cet exemple d'erreur](https://stackoverflow.com/questions/26920083/fatal-error-mpi-h-no-such-file-or-directory-include-mpi-h))
- Positionner la variable `HYPRE_DIR` dans le fichier `Makefile`
- Lancer la compilation avec `make tp`  
