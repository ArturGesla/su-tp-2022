#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 4
#define EPS 1.e-8

double L_oo(double x[])
{
  double m = fabs(x[0]);
  for (int i = 1; i < N; i++)
    if (fabs(x[i]) > m)
      m = fabs(x[i]);
  return m;
}

double prod_scal(double x[], double y[])
{
  double x_scal_y = 0;
  for (int i = 0; i < N; i++)
    x_scal_y += x[i] * y[i];
  return x_scal_y;
}

void prod_Ax(double A[N][N], double x[], double y[])
{
  for (int i = 0; i < N; i++)
  {
    y[i] = 0;
    for (int j = 0; j < N; j++)
    {
      y[i] += A[i][j] * x[j];
    }
  }
}
void prod_cx(double c, double x[], double y[])
{
  for (int i = 0; i < N; i++)
  {
    y[i] = c * x[i];
  }
}
double L_2(double x[])
{
  return sqrt(prod_scal(x, x));
}

void transpose(double A[][N], double At[][N])
{
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      At[i][j] = A[j][i];
  return;
}

void afficher_vect(double x[])
{
  for (int i = 0; i < N; i++)
    printf("%1lf\n", x[i]);
  return;
}

void afficher_mat(double A[][N])
{
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
      printf("%.4lf\t", A[i][j]);
    printf("\n");
  }
  return;
}

// DESCENTE
void resol_trig_inf(double L[][N], double x[], double b[])
{
  for (int i = 0; i < N; i++)
  {
    double l = 0;
    for (int k = 0; k < i; k++)
      l += L[i][k] * x[k];
    x[i] = (b[i] - l) / L[i][i];
  }
  return;
}

// MONTÉE
void resol_trig_sup(double U[][N], double x[], double b[])
{
  for (int i = N - 1; i >= 0; i--)
  {
    double l = 0;
    for (int k = i + 1; k < N; k++)
      l += U[i][k] * x[k];
    x[i] = (b[i] - l) / U[i][i];
  }
  return;
}

void fact_LU(double A[][N], double L[][N], double U[][N])
{
  for (int j = 0; j < N; j++)
  {
    for (int i = 0; i <= j; i++)
    {
      double l = 0;
      for (int k = 0; k < i; k++)
        l += L[i][k] * U[k][j];
      U[i][j] = A[i][j] - l;
    }
    L[j][j] = 1;
    for (int i = j + 1; i < N; i++)
    {
      double l = 0;
      for (int k = 0; k < j; k++)
        l += L[i][k] * U[k][j];
      L[i][j] = (A[i][j] - l) / U[j][j];
    }
  }
  return;
}

double puiss_it(double A[][N], double x[])
{
  double u[N], lambda, lambda_old;
  double y[N];
  lambda = 1;
  // afficher_vect(y);
  // A FAIRE!

  int nit = 40;
  // for (int i = 0; i < nit; i++)
  int i = 0;
  while (fabs(lambda_old - lambda) > 1e-15 && i < nit)
  {
    prod_Ax(A, x, y);
    lambda_old = lambda;
    lambda = L_2(y);
    // printf("it: %d current ev: %4.8f dev: %4.2e \n", i, lambda, fabs(lambda_old - lambda));
    // normalise
    prod_cx(1 / lambda, y, x);
    // printf("current evc: \n");
    // afficher_vect(x);
    i++;
  }
  printf("it: %d current ev: %4.8f dev: %4.2e \n", i, lambda, fabs(lambda_old - lambda));

  if (abs(lambda_old - lambda) > 1e-16)
  {
    printf("========== Error: no convergence. ==========");
  }
  return lambda;
}

double puiss_inv(double A[][N], double x[])
{
  double L[N][N], U[N][N];
  double v[N], y[N];
  double lambda, lambda_old;
  lambda=1;

  int nit = 40;
  // for (int i = 0; i < nit; i++)
  int i = 0;
      fact_LU(A, L, U);

  while (fabs(lambda_old - lambda) > 1e-15 && i < nit)
  {
    resol_trig_inf(L, v, x);
    resol_trig_sup(U, y, v);
    lambda_old = lambda;
    lambda = L_2(y);
    printf("it: %d current ev: %4.8f dev: %4.2e \n", i, lambda, fabs(lambda_old - lambda));
    prod_cx(1 / lambda, y, x);
    i++;
  }
  printf("it: %d current ev: %4.8f dev: %4.2e \n", i, lambda, fabs(lambda_old - lambda));

  if (abs(lambda_old - lambda) > 1e-16)
  {
    printf("========== Error: no convergence. ==========");
  }
  return lambda;

  // A FAIRE!
}

void deflation_it(double A[][N], double x[], double lambda)
{
  double At[N][N], xt[N];
  prod_cx(1, x, xt);
  transpose(A, At);
  double l1 = puiss_it(At, xt);
  if (fabs(l1 - lambda) > 1e-14)
    printf("========== Error: diff ev of transpose. Diff: %4.2e==========\n", fabs(l1 - lambda));

  // construct deflated operator
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      A[i][j] = A[i][j] - lambda * x[i] * xt[j] / prod_scal(x, xt);
    }
  }
  // afficher_mat(A);
  // afficher_mat(At);

  // A FAIRE!
}

int main(void)
{
  double x1[N], x2[N], x3[N], x4[N];
  double A[N][N];
  double L1, L2, L3, L4;

  // Initialisation de A ici matrice de Vandermonde - A FAIRE

  // diffusion
  {
    int n = 3;
    for (int i = 0; i < n; i++)
    {
      A[i][i] = 2.0;
      if (i > 0)
        A[i][i - 1] = -1;
      if (i < n - 1)
        A[i][i + 1] = -1;
    }

    afficher_mat(A);
    double x1[N] = {0, 0, 1};
    prod_cx(1, x1, x2);
    prod_cx(1, x1, x3);

    // double x0[N] = {1.0 / 2, sqrt(2) / 2, 1.0 / 2};
    // afficher_vect(x3);

    // afficher_vect(x1);
    L1 = puiss_it(A, x1);
    deflation_it(A, x1, L1);
    L2 = puiss_it(A, x2);
    deflation_it(A, x2, L2);
    L3 = puiss_it(A, x3);
  }

  // Vandermonde
  {
    int n = 4;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        A[i][j] = pow((i + 1), j + 1);
      }
    }
    afficher_mat(A);
    double x1[N] = {0, 0, 1};
    prod_cx(1, x1, x2);
    prod_cx(1, x1, x3);
    prod_cx(1, x1, x4);

    // double x0[N] = {1.0 / 2, sqrt(2) / 2, 1.0 / 2};
    // afficher_vect(x3);

    // afficher_vect(x1);

    // L1 = puiss_it(A, x1);
    // deflation_it(A, x1, L1);
    // L2 = puiss_it(A, x2);
    // deflation_it(A, x2, L2);
    // L3 = puiss_it(A, x3);
    // deflation_it(A, x3, L3);
    // L4 = puiss_it(A, x4);

    L1 = puiss_inv(A, x1);
  }

  // Recherche des valeurs propres par la puissance itérée

  return 0;
}
