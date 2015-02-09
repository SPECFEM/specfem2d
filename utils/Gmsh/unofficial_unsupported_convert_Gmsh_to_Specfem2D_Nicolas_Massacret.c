#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define TAILLE_MAX 100 // Tableau de taille 100

int main(int argc, char *argv[])
{
  FILE* fichier = NULL;
      FILE* fichierNode = NULL;
  FILE* fichierElem = NULL;
  FILE* fichierMat = NULL;
      char chaine[TAILLE_MAX] = "";
  char nodes[TAILLE_MAX] = "$Nodes";
      char endnodes[TAILLE_MAX] = "$EndNodes";
  char elem[TAILLE_MAX] = "$Elements";
  float coordNoeud[3] = {0};
  int car = 0;
  int x=0;
  long infos[2] = {0};
  long numNode[4] = {0};
  long i;
  char nomFichier[20] = "";


// printf("-> Indiquez le nom du fichier source : ");
// scanf("%s", &nomFichier[0]);
// printf("\n - Ouvertur du fichier source du maillage : %s \n", nomFichier[0]);

      fichier = fopen("jpz1.msh", "r");
  if (fichier != NULL)
  {
    fscanf(fichier, "%s", &chaine[0]);
    while ( strcmp(&chaine[0], &nodes[0]) != 0)
    {
      fscanf(fichier, "%s", &chaine[0]);
    }             // Lorsque l'on sort de cette boucle, on a trouve $Nodes

  /* FDEBUT DE L'ECRITURE DE NODE */
    fichierNode = fopen("nodes", "r+");
    if (fichierNode != NULL)
        {
      fscanf(fichier, "%li", &infos[0]);
      fprintf(fichierNode, "%li \n", infos[0]);     // on note le nombre de noeuds
int nbrNode = infos[0];
      car = fgetc(fichier);
      while(  car != 36)        // Boucle permettant l'ecriture du fichier contenant les coordonees des noeuds?.
      {
        fseek(fichier, -1, SEEK_CUR);
        fscanf(fichier, "%f %f %f %f", &coordNoeud[0], &coordNoeud[1], &coordNoeud[2], &coordNoeud[3]);
        fprintf(fichierNode, "%f %f \n", coordNoeud[1], coordNoeud[2]);

      if(coordNoeud[2] == 0)    // Calcul du nombre d'elements suivant x.
      {
        x++;
      }

int etat = ((coordNoeud[0]*100) / nbrNode);
printf("\r- Ecriture des coordonees des noeuds : %d %%", etat);
        car = fgetc(fichier);
        car = fgetc(fichier);
      }

      fclose(fichierNode);
    }
    else {printf("!!! Erreur : Probleme a l'ouverture du fichier nodes. \n"); return 0;}
 /* FIN DE L'ECRITURE DE NODE */

    printf("\n- Ecriture du fichier node terminee. \n");
    x--;
    printf("- Nombre d'elements suivant x = %d \n", x);

/* DEBUT DE L'ECRITURE DE ELEMENT ET METERIAL */
    fichierElem = fopen("element", "r+");
    if (fichierElem != NULL)
          {
// printf("debut ecriture element et material \n");
      fichierMat = fopen("material", "r+");
      if (fichierMat != NULL)
        {
        fscanf(fichier, "%s", &chaine[0]); // pour sauter une ligne
        fscanf(fichier, "%s", &chaine[0]);
        fscanf(fichier, "%li", &infos[0]);
        fprintf(fichierElem, "%li \n", infos[0]); // on note le nombre d'elements
int nbrElem = infos[0];
        car = fgetc(fichier);
        while(  car != 36)      // Pour s'arreter au prochain $.
          {
          fseek(fichier, -1, SEEK_CUR);
          fscanf(fichier, "%li %li %li", &infos[0], &infos[1], &infos[2]);
          if (infos[1] == 3 )   // Pour ne retenir que les maillage hexa.
            {
            int j = infos[2];
            for(i=0; i<j; i++) //on saute les tags
              {
              fscanf(fichier, "%li", &numNode[i]);
              }
int etat = ((infos[0]*100) / nbrElem);
printf("\r- Ecriture des elements et materiaux : %d %%", etat);
            fprintf(fichierMat, "%li \n", numNode[0]);    //on note le numero du matos
            fscanf(fichier, "%li %li %li %li", &numNode[0], &numNode[1], &numNode[2], &numNode[3]);
            fprintf(fichierElem, "%li %li %li %li\n", numNode[0], numNode[1], numNode[2], numNode[3]); // on note les ID des nodes
            }

          car = fgetc(fichier);
          car = fgetc(fichier);
          }

        fclose(fichierMat);
        }

      else
        {
        printf("Erreur ouverture Material \n"); return 0;
        }

      fclose(fichierElem);
      }

    else
      {
      printf("Erreur ouverture Element \n"); return 0;
      }

/* FIN DE L'ECRITURE DE ELEMENT ET MATERIAL */

        fclose(fichier);
      }
else
  {
  printf("Erreur ouverture fichier \n"); return 0;
  }
printf("\n- Ecriture des fichiers element et material terminee. \n");
return 0;
}

