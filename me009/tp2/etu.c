// Partie II - Conditions Mixtes
// Natacha Guegan-Fau 28709310
// Anna Sizonenko 28707738
// Groupe 6A

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#define MAX 35000
#define NMAX 100

void afficher_vect(double x[MAX], int n)
{
    int i;
    for (i = 0; i < n; i++)
        printf("%d\t%lf\n", i, x[i]);
}

double norme2_vect(double x[MAX], int n)
{
    double var = 0;
    int i;
    for (i = 0; i < n; i++)
        var += (x[i] * x[i]);
    return sqrt(var);
}

void ecriture(double Q[MAX], double P[NMAX][NMAX], double x[NMAX], double y[NMAX], int nx, int ny, double h, double gw)
{
    int i, j, JJ;

    FILE *fichier1;
    fichier1 = fopen("isoP_sor.dat", "w");

    // remplissage des coordonnées du maillage x, y et de solution P en chaque point du maillage
    for (i = 0; i <= nx + 1; i++)
    {
        x[i] = i * h;
    }
    for (i = 0; i <= ny + 1; i++)
    {
        y[i] = i * h;
    }
    puts("Remplissage du maillage REUSSIT");

    for (j = 0; j <= ny; j++)
    {
        P[0][j] = gw;
    }

    // points aux sommets
    P[0][0] = 1. / 0.;
    P[0][ny + 1] = 1. / 0.;
    P[nx + 1][0] = 1. / 0.;
    P[nx + 1][ny + 1] = 1. / 0.;
    puts("Application des conditions limites REUSSIT");

    for (i = 1; i <= nx; i++)
    {
        for (j = 0; j <= ny + 1; j++)
        {
            JJ = i + j * nx - 1;
            P[i][j] = Q[JJ];
        }
    }

    for (j = 0; j <= ny + 1; j++)
    {
        P[nx + 1][j] = y[j];
    }

    if (nx == 5)
    {
        for (i = 0; i <= nx + 1; i++)
            printf("%e ", x[i]);
        printf("\n");
        for (i = 0; i <= ny + 1; i++)
            printf("%e ", y[i]);
        printf("\n");
        puts("P");

        for (j = ny + 1; j >= 0; j--)
        {
            for (i = 0; i <= nx + 1; i++)
            {
                printf("%e ", P[i][j]);
            }
            printf("\n");
        }
    }

    // Ecriture du maillage et du résultat dans un fichier
    for (j = 0; j <= ny + 1; j++)
    {
        for (i = 0; i <= nx + 1; i++)
        {
            fprintf(fichier1, "%lf\t%lf\t%lf\n", x[i], y[j], P[i][j]);
        }
        fprintf(fichier1, "\n");
    }
    fclose(fichier1);
    puts("Creation du fichier avec maillage et resultats REUSSIT\n");
}

int main()
{
    int nx, ny, nn;
    int c;
    printf("Choisissez le maillage:\n1)  nx=5  ny=4\n2)  nx=47  ny=39\n3)  nx=95  ny=79\n4)  nx=191  ny=159\n");
    scanf("%d", &c);
    if (c == 1)
    {
        nx = 5, ny = 4, nn = nx * ny;
    }
    else if (c == 2)
    {
        nx = 47, ny = 39, nn = nx * ny;
    }
    else if (c == 3)
    {
        nx = 95, ny = 79, nn = nx * ny;
    }
    else if (c == 4)
    {
        nx = 191, ny = 159, nn = nx * ny;
    }
    else
    {
        printf("Erreur de saisie\n");
    }
    printf("Maillage choisi: nx = %d, ny = %d\n", nx, ny);

    double An[MAX], Ap[MAX], As[MAX], Ae[MAX], Aw[MAX], b[MAX], r[MAX];
    double h = 1.5 / (nx + 1);
    double gw = 0., f = 0;
    double ge = 0;

    // Pour la résolution SOR
    int kmax = 5000, kech = 20;
    double eps = 1.e-8;
    double omega, omega_opt = 0., omega_min = 1., omega_max = 2.;
    double nores, nores_min = 1.e8;
    double Q[MAX], dQ[MAX];

    // Champ solution 2D
    double P[NMAX][NMAX];

    // Coordonnées du maillage
    double x[NMAX], y[NMAX];

    // Pour le calcul de flux
    int mx = (nx + 1) / 2, my = (ny + 1) / 2;
    double T1 = 293., T2 = 333, lambda = 120;
    double gtv[NMAX], gth[NMAX];
    double IH, IV, QH, QV;

    int i, j, k;

    // implémentation des vecteurs An,Ap,As,Ae,Aw
    for (int i = 0; i < nx; i++)
    {
        Ap[i] = 2;
        Aw[i] = Ae[i] = -1. / 2;
        An[i] = -1;
        As[i] = 0;
        b[i] = h * h * f;
    }
    Aw[0] = Ae[nx - 1] = 0;

    for (int i = 0 + nx; i < nn + nx; i++)
    {
        Ap[i] = 4;
        An[i] = Aw[i] = As[i] = Ae[i] = -1;
        b[i] = h * h * f;
    }

    for (int i = nx; i < nn + nx; i = i + nx)
    {
        b[i] = b[i] + gw + h * h * f;
        Aw[i] = 0;
    }

    for (int i = 2 * nx - 1; i < nn + nx; i = i + nx)
    {
        ge = ge + h;
        b[i] = b[i] + ge + h * h * f;
        Ae[i] = 0;
    }

    for (int i = nn + nx; i < nn + 2 * nx; i++)
    {
        Ap[i] = 2;
        Aw[i] = Ae[i] = -1. / 2;
        As[i] = -1;
        An[i] = 0;
        b[i] = h * h * f;
    }

    Aw[nn + nx] = Ae[nn + 2 * nx - 1] = 0;
    ge = ge + h;
    // ge = ge + h/2.0; // maybe
    b[nn + 2 * nx - 1] = ge + h * h * f;
    puts("\nImplémentation des vecteurs An,Ap,As,Ae,Aw et b REUSSITE");

    FILE *fichier = NULL;
    fichier = fopen("res_sor.dat", "w");

    // résolution par SOR
    puts("\nMéthode SOR");

    // choix de omega optimal
    // for (omega = omega_min; omega <= omega_max; omega += 0.01)
        omega =1.91;

    {
        for (int i = 0; i < nn; i++)
        {
            Q[i] = 0.;
        }
        for (int i = 0; i < nn; i++)
        {
            r[i] = b[i];
        }
        k = 0;

        while (norme2_vect(r, nn) > eps && k < kmax)
        {
            r[0] = b[0] - Ap[0] * Q[0] - Ae[0] * Q[1] - An[0] * Q[nx];
            dQ[0] = omega * r[0] / Ap[0];

            for (i = 1; i < nx; i++)
            {
                r[i] = b[i] - Aw[i] * Q[i - 1] - Ap[i] * Q[i] - Ae[i] * Q[i + 1] - An[i] * Q[i + nx];
                dQ[i] = omega * (r[i] - Aw[i] * dQ[i - 1]) / Ap[i];
            }

            for (i = nx; i < nn + nx; i++)
            {
                r[i] = b[i] - As[i] * Q[i - nx] - Aw[i] * Q[i - 1] - Ap[i] * Q[i] - Ae[i] * Q[i + 1] - An[i] * Q[i + nx];
                dQ[i] = omega * (r[i] - As[i] * dQ[i - nx] - Aw[i] * dQ[i - 1]) / Ap[i];
            }

            for (i = nn + nx; i < nn + 2 * nx - 1; i++)
            {
                r[i] = b[i] - As[i] * Q[i - nx] - Aw[i] * Q[i - 1] - Ap[i] * Q[i] - Ae[i] * Q[i + 1];
                dQ[i] = omega * (r[i] - As[i] * dQ[i - nx] - Aw[i] * dQ[i - 1]) / Ap[i];
            }

            r[nn + 2 * nx - 1] = b[nn + 2 * nx - 1] - As[nn + 2 * nx - 1] * Q[nn + nx - 1] - Aw[nn + 2 * nx - 1] * Q[nn + 2 * nx - 2] - Ap[i] * Q[i];
            dQ[nn + 2 * nx - 1] = omega * (r[nn + 2 * nx - 1] - As[nn + 2 * nx - 1] * dQ[nn + nx - 1] - Aw[i] * dQ[nn + 2 * nx - 2]) / Ap[nn + 2 * nx - 1];

            for (i = 0; i < nn + 2 * nx; i++)
            {
                Q[i] = Q[i] + dQ[i];
            };

            nores = norme2_vect(r, nn);
            k = k + 1;
        }

        // mis à jour du facteur omega
        if (nores < nores_min)
        {
            nores_min = nores;
            omega_opt = omega;
        }
    }

    printf("omega_opt = %lf\n\n", omega_opt);

    // resolution avec omega optimal
    puts("méthode SOR");
    for (int i = 0; i < nn + 2 * nx; i++)
    {
        Q[i] = 0;
    }
    for (int i = 0; i < nn + 2 * nx; i++)
    {
        r[i] = b[i];
    }
    k = 0;

    while (norme2_vect(r, nn) > eps && k < kmax)
    {

        r[0] = b[0] - Ap[0] * Q[0] - Ae[0] * Q[1] - An[0] * Q[nx];
        dQ[0] = omega_opt * r[0] / Ap[0];

        for (i = 1; i < nx; i++)
        {
            r[i] = b[i] - Aw[i] * Q[i - 1] - Ap[i] * Q[i] - Ae[i] * Q[i + 1] - An[i] * Q[i + nx];
            dQ[i] = omega_opt * (r[i] - Aw[i] * dQ[i - 1]) / Ap[i];
        }

        for (i = nx; i < nn + nx; i++)
        {
            r[i] = b[i] - As[i] * Q[i - nx] - Aw[i] * Q[i - 1] - Ap[i] * Q[i] - Ae[i] * Q[i + 1] - An[i] * Q[i + nx];
            dQ[i] = omega_opt * (r[i] - As[i] * dQ[i - nx] - Aw[i] * dQ[i - 1]) / Ap[i];
        }

        for (i = nn + nx; i < nn + 2 * nx - 1; i++)
        {
            r[i] = b[i] - As[i] * Q[i - nx] - Aw[i] * Q[i - 1] - Ap[i] * Q[i] - Ae[i] * Q[i + 1];
            dQ[i] = omega_opt * (r[i] - As[i] * dQ[i - nx] - Aw[i] * dQ[i - 1]) / Ap[i];
        }

        r[nn + 2 * nx - 1] = b[nn + 2 * nx - 1] - As[nn + 2 * nx - 1] * Q[nn + nx - 1] - Aw[nn + 2 * nx - 1] * Q[nn + 2 * nx - 2] - Ap[i] * Q[i];
        dQ[nn + 2 * nx - 1] = omega_opt * (r[nn + 2 * nx - 1] - As[nn + 2 * nx - 1] * dQ[nn + nx - 1] - Aw[i] * dQ[nn + 2 * nx - 2]) / Ap[nn + 2 * nx - 1];

        for (i = 0; i < nn + 2 * nx; i++)
        {
            Q[i] = Q[i] + dQ[i];
        }
        nores = norme2_vect(r, nn);

        if ((k % kech) == 0)
        {
            // printf("k = %d, résidu = %e\n", k, nores);
            fprintf(fichier, "%d\t%e\n", k, nores);
            // k désigne ici le nombre d'itérations, nores désigne ici la norme du résidue
        }
        k = k + 1;
    }

    printf("converegnce en k = %d iterations, résidu = %e\n", k, nores);
    fprintf(fichier, "%d\t%e\n", k, nores);

    puts("");
    puts("vecteur Q pour nx==5");
    if (nx == 5)
    {
        for (i = 0; i < nn + 2 * nx; i++)
            printf(" %e", Q[i]);
        printf("\n");
    }

    puts("écriture");
    ecriture(Q, P, x, y, nx, ny, h, gw);
    // Calcul de flux QH et QV
    IH = 0.;
    IV = 0.;

    for (i = 0; i <= nx; i++)
    {
        gtv[i] = -P[i][my + 1] - P[i + 1][my + 1] + P[i][my - 1] + P[i + 1][my - 1];
        IV = IV + gtv[i];
    }
    IV = 1 / 4. * IV;

    for (i = 0; i <= ny; i++)
    {
        gth[i] = -P[mx + 1][i] + P[mx + 1][i + 1] - P[mx - 1][i] + P[mx - 1][i + 1];
        IH = IH + gth[i];
    }
    IH = 1 / 4. * IH;

    QV = lambda * (T2 - T1) * IV;
    QH = lambda * (T2 - T1) * IH;

    printf("IH = %lf, IV = %lf\n", IH, IV);
    printf("QH = %lf, QV = %lf\n", QH, QV);
    fclose(fichier);

    // Solution exacte
    puts("\nCalculs des ecarts");
    if (nx == 95)
    {
        double T[NMAX][NMAX] = {0};
        double lambda, A, err_global, err, somme_err;

        int k, l;
        char nom_fichier[19] = "isoP_exacte_";
        char fichier_err[19] = "erreur_N_";
        FILE *fichier;
        FILE *fichier2;
        int Nserie[4] = {10, 25, 50, 100};
        for (k = 0; k < 4; k++)
        {

            sprintf(&nom_fichier[12], "%d.dat", Nserie[k] * 2);
            sprintf(&fichier_err[9], "%d.dat", Nserie[k] * 2);
            fichier = fopen(nom_fichier, "w+");
            fichier2 = fopen(fichier_err, "w+");
            if (errno)
            {
                printf("%s\n", strerror(errno));
                return (1);
            }
            // calculs solution approche, erreur et sauvegarde dans fihiers
            somme_err = 0;
            for (i = 0; i <= nx + 1; i++)
            {
                for (j = 0; j <= ny + 1; j++)
                {
                    T[i][j] = x[i] * 1.25 / (2 * 1.5);
                    for (l = 0; l < Nserie[k]; l++)
                    {
                        lambda = (2 * l + 1) * M_PI / 1.25;
                        A = -4 * 1.25 / ((2 * l + 1) * (2 * l + 1) * M_PI * M_PI * sinh(lambda * 1.5));
                        T[i][j] = T[i][j] + A * sinh(lambda * x[i]) * cos(lambda * y[j]);
                    }
                    fprintf(fichier, "%lf\t%lf\t%lf\n", x[i], y[j], T[i][j]);
                    err = fabs(T[i][j] - P[i][j]);
                    fprintf(fichier2, "%lf\t%lf\t%lf\n", x[i], y[j], err);
                    if ((i != 0) & (j != 0))
                    {
                        somme_err = somme_err + err * err;
                    }
                }
            }
            fclose(fichier);
            fclose(fichier2);

            // calcul de l'erreur global
            err_global = sqrt(somme_err / (nx * ny));
            printf("N = %d\tsol_app = %lf\terr = %lf\terr global = %lf\n\n", Nserie[k] * 2, T[nx + 1][ny + 1], err, err_global);
        }
    }
    return 0;
}