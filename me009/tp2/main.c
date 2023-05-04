#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"
#define MAX 160000
#define NMAX 400

float norme2_vect(float x[MAX], int n)
{
    float var = 0;
    int i;
    for (i = 0; i < n; i++)
        // var += (x[i] * x[i]);
        var += (x[i] * x[i] / n);

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
            fprintf(fichier1, "%4.8f\t", Q[nx * iy + ix]);
        }
        fprintf(fichier1, "\n");
    }

    fclose(fichier1);
}

int main()
{

    // Discretisation

    // int nx = 5, ny = 4, nn = nx * ny;

    // int nx = 11, ny = 9, nn = nx * ny;

    // int nx = 23, ny = 19, nn = nx * ny;

    // int nx = 47, ny = 39, nn = nx * ny;

       int nx = 95, ny = 79, nn = nx * ny;
    //    three above give ration 3.72810099409025 of QV, not too bad

    //  int nx = 191, ny = 159, nn = nx * ny;
    //   int nx = 383, ny = 319, nn = nx * ny;

    float h = 1.5 / (nx + 1);

    // Conditions aux limites
    float gs = 0., gw = 1., ge = 0.2, gn = 0.6, f = 0;
    // float gs = 0.3, gw = 0.5, ge = 1., gn = 0.0;

    // Diagonales non nulles de la matrice, vecteur second membre et r�sidu
    float As[MAX], Aw[MAX], Ap[MAX], Ae[MAX], An[MAX], b[MAX], r[MAX];

    // Pour la r�solution SOR
    int kmax = 15000, kech = 20;
    float eps = 2e-7, omega = 1.8, nores = 1.;
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
    // remplirDirichlet(An, As, Ap, Ae, Aw, b, nn, nx, ny, gs, ge, gw, gn);
    ny = ny + 2;
    nn = nx * ny;
    remplirNeumann(An, As, Ap, Ae, Aw, b, nn, nx, ny, gs, ge, gw, gn, h);

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
            printf("%4.4f, ", b[i]);
        printf("\n");
    }

    // R�solution par SOR

    // Initialisation
    initQ(Q, nn);
    // A compl�ter

    // Boucle conditionnelle
    while (nores > eps /*� completer*/ && k < kmax)
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
            fprintf(fichier, "%d\t %e\n", k, nores);
            /*k d�signe ici le nombre d'it�rations, nores d�signe ici la norme du r�sidue*/
        }

        // A compl�ter
        k++;
    }
    fprintf(fichier, "%d\t %e\n", k, nores);

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
    printf("QH = %f, QV = %2.8e\n", QH, fluxV4Dirichlet(Q, nx, ny, h));

    fclose(fichier);
    int nmodes = 60;
    float nmDiff=diffWithAnalytical(Q, nn, nx, ny, h, nmodes);
    ecriture(Q, P, x, y, nx, ny, h, gs, gw, ge, gn);
    printf("norm of diff with analytical %4.2e\n", nmDiff);

    return 0;
}
