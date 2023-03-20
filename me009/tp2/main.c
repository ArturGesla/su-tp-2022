#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#define MAX 10000
#define NMAX 100

float norme2_vect(float x[MAX], int n)
{
    float var = 0;
    int i;
    for (i = 0; i < n; i++)
        var += (x[i] * x[i]);

    return sqrt(var);
}

void ecriture(float Q[MAX], float P[NMAX][NMAX], float x[NMAX], float y[NMAX],
              int nx, int ny, float h, float gs, float gw, float ge, float gn)
{

    int i, j, JJ;

    FILE *fichier1;
    fichier1 = fopen("isoP_sor.dat", "w");

    // Remplissage des coordonn�es du maillage x, y et de solution P en chaque point du maillage

    // A compl�ter

    if (nx == 5)
    {
        for (i = 0; i <= nx + 1; i++)
            printf("%e ", x[i]);
        printf("\n");

        for (i = 0; i <= ny + 1; i++)
            printf("%e ", y[i]);
        printf("\n");

        for (i = 0; i <= nx + 1; i++)
        {
            for (j = 0; j <= ny + 1; j++)
            {
                printf("%e ", P[i][j]);
            }
            printf("\n");
        }
    }

    // Ecriture du maillage et du r�sultat dans un fichier

    // � compl�ter
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iy = 0; iy < ny; iy++)
        {
            fprintf(fichier1,"%4.8f\t",Q[nx*iy+ix]);
        }
        fprintf(fichier1,"\n");
        
    }
    

    fclose(fichier1);
}

int main()
{

    // Discretisation
    int nx = 50, ny = 40, nn = nx * ny;
    //  int nx = 30, ny = 25, nn = nx * ny;
    //    int nx = 95, ny = 79, nn = nx * ny;
    float h = 1.5 / (nx + 1);

    // Conditions aux limites
    float gs = 0., gw = 1., ge = 0.2, gn = 0.6, f = 0;
    // float gs = 0.3, gw = 0.5, ge = 1., gn = 0.0;

    // Diagonales non nulles de la matrice, vecteur second membre et r�sidu
    float As[MAX], Aw[MAX], Ap[MAX], Ae[MAX], An[MAX], b[MAX], r[MAX];

    // Pour la r�solution SOR
    int kmax = 5000, kech = 20;
    float eps = 1.e-6, omega = 1.9, nores = 1.;
    float Q[MAX], dQ[MAX];

    // Champ solution 2D
    float P[NMAX][NMAX];

    // Coordonn�es du maillage
    float x[NMAX], y[NMAX];

    // Pour le calcul de flux
    int mx, my;
    float T1 = 293., T2 = 333., lambda = 120;
    float gtv[NMAX], gth[NMAX];
    float IH, IV, QH, QV;

    int i, j, k, JJ;

    // Mise en equation
    // initiation des vecteurs non nuls de la matrice pentadiagonale As, Aw, Ap, Ae, An
    // et du second membre

    // A compl�ter
    remplirDirichlet(An, As, Ap, Ae, Aw, b, nn, nx, ny, gs, ge, gw, gn);
    // Prise en compte des conditions aux limites

    // A compl�ter

    FILE *fichier = NULL;
    fichier = fopen("res_sor.dat", "w");

    printf("nx = %d, ny = %d\n", nx, ny);

    // impression de la matrice et du second membre
    if (nx == 5)
    {
        for (int i = 0; i < nn; i++)
            printf("%.1f, ", As[i]);
        printf("\n");
        for (int i = 0; i < nn; i++)
            printf("%.1f, ", Aw[i]);
        printf("\n");
        for (int i = 0; i < nn; i++)
            printf("%.1f, ", Ap[i]);
        printf("\n");
        for (int i = 0; i < nn; i++)
            printf("%.1f, ", Ae[i]);
        printf("\n");
        for (int i = 0; i < nn; i++)
            printf("%.1f, ", An[i]);
        printf("\n");
        for (int i = 0; i < nn; i++)
            printf("%.1f, ", b[i]);
        printf("\n");
    }

    // R�solution par SOR

    // Initialisation
    initQ( Q,nn);
    // A compl�ter

    // Boucle conditionnelle
    while (nores > eps /*� completer*/ && k<5000)
    {

        // A compl�ter
        bMinusAq(An, As, Ap, Ae, Aw, b, nn, nx, ny, Q, r);
        descentedq(An, As, Ap, Ae, Aw, b, nn, nx, ny, r, dQ, omega);
        amendQ(dQ, Q, nn);
        nores = norme2_vect(r, nn);

        if ((k % kech) == 0)
        {
            printf("k = %d, r�sidu = %e\n", k, nores);
            // fprintf(fichier, "k = %d, r�sidu = %e\n", k, nores);
            fprintf(fichier,"%d\t %e\n", k, nores);
            /*k d�signe ici le nombre d'it�rations, nores d�signe ici la norme du r�sidue*/
        }

        // A compl�ter
        k++;
    }
            fprintf(fichier,"%d\t %e\n", k, nores);

    printf("converegnce en k = %d iterations, r�sidu = %e\n", k, nores);
    // fprintf(fichier, "converegnce en k = %d iterations, r�sidu = %e\n", k, nores);

    if (nx == 5)
    {
        for (i = 0; i < nn; i++)
            printf(" %e", Q[i]);
        printf("\n");
    }

    ecriture(Q, P, x, y, nx, ny, h, gs, gw, ge, gn);

    // Calcul de flux QH et QV

    // A compl�ter

    printf("IH = %f, IV = %f\n", IH, IV);
    printf("QH = %f, QV = %f\n", QH, QV);

    fclose(fichier);

    return 0;
}
