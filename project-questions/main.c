#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 500
#include <time.h>

// tp1
void factoriser_LU(const float A[MAX][MAX], float L[MAX][MAX], float U[MAX][MAX], int n)
{
	//	int i,j,k;
	// factorisation à faire
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < j + 1; i++)
		{
			U[i][j] = A[i][j];
			// printf("%4.2f \t ",U[i][j]);
			// printf("%4.2f \n ",A[i][j]);
			for (int k = 0; k < i; k++)
			{
				U[i][j] -= L[i][k] * U[k][j];
			} // k
		}	  // i
		L[j][j] = 1;
		for (int i = j + 1; i < n; i++)
		{
			L[i][j] = A[i][j];
			// printf("%4.2f \n ",L[i][j]);

			for (int k = 0; k < j; k++)
			{
				L[i][j] -= L[i][k] * U[k][j];
			} // k
			L[i][j] = L[i][j] / U[j][j];
		} // i

	} // j
}

void resol_trig_inf(const float A[MAX][MAX], float x[MAX], const float b[MAX], int n)
{
	//	int i,k;

	// la descente à faire
	for (int i = 0; i < n; i++)
	{
		x[i] = b[i];
		for (int j = 0; j < i; j++)
		{
			x[i] -= x[j] * A[i][j];
		}
		x[i] = x[i] / A[i][i];
	}
}

void resol_trig_sup(float A[MAX][MAX], float x[MAX], float b[MAX], int n)
{
	// la remontée à faire
	for (int i = n - 1; i > -1; i--)
	{
		x[i] = b[i];
		for (int j = n - 1; j > i; j--)
		{
			x[i] -= x[j] * A[i][j];
		}
		x[i] = x[i] / A[i][i];
	}
}

void afficher_vect(float x[MAX], int n)
{
	int i;
	for (i = 0; i < n; i++)
		printf("%4.4f\n", x[i]);
}

void afficher_mat(float A[MAX][MAX], int n)
{
	int i, j;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%.1f\t", A[i][j]);
		}
		printf("\n");
	}
}

void remplir_vect(float x[MAX], int n)
{
	int i;
	for (i = 0; i < n; i++)
	{
		printf("valeur [%d] : ", i + 1);
		if (scanf("%f", &x[i]) != 1)
			printf("error");
	}
}

void remplir_mat(float A[MAX][MAX], int n)
{
	int i, j;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("valeur [%d][%d] : ", i + 1, j + 1);
			if (scanf("%f", &A[i][j]) != 1)
				printf("error");
		}
		printf("\n");
	}
}

void remplir_mat_rand(float A[MAX][MAX], int n)
{

	srand(time(NULL));
	int i, j;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			float r = rand();
			A[i][j] = r / RAND_MAX;
		}
	}
}

// tp2
void calculateResid(const float A[MAX][MAX], const float b[], const float u[], float r[], int n)
{
	for (int i = 0; i < n; i++)
	{
		r[i] = 0;
		for (int j = 0; j < n; j++)
		{
			r[i] += -A[i][j] * u[j];
		}
		r[i] += b[i];
	}
}

void calculateM(const float A[MAX][MAX], float M[MAX][MAX], int n)
{
	float omega = 1.0;
	// float omega = 1.0;

	for (int i = 0; i < n; i++)
	{
		// M[i][i] = 1; //naive
		M[i][i] = A[i][i] / omega; // jacobi
	}

	// gs
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			//  M[i][j] = +A[i][j];
		}
	}
}

void calculateDu(const float M[MAX][MAX], const float r[MAX], float du[MAX], int n)
{
	resol_trig_inf(M, du, r, n);
}

void incrementU(const float du[MAX], float u[MAX], int n)
{
	for (int i = 0; i < n; i++)
	{
		u[i] = u[i] + du[i];
	}
}
float norm(const float r[MAX], int n)
{
	float result = 0;
	for (int i = 0; i < n; i++)
	{
		result += r[i] * r[i];
	}

	return sqrtf(result);
}

void solveIter(const float A[MAX][MAX],
			   const float b[MAX], float u[MAX],
			   int n, float eps, int nsteps)
{
	float r[MAX];
	float M[MAX][MAX];
	calculateM(A, M, n);
	float resN;
	// float resNm1;
	// resN=norm(r, n);
	float normDu = 1;
	float normDum1;

	for (int isteps = 0; isteps < nsteps; isteps++)
	{
		calculateResid(A, b, u, r, n);
		// resNm1=resN;
		resN = norm(r, n);
		// printf("iter: %d \t resid: %4.2e \n", isteps, resN/resNm1);
		printf("iter: %d \t resid: %4.2e \t", isteps, resN);

		if (resN < 1e-4)
			break;
		// afficher_vect(r, n);
		float du[MAX];
		calculateDu(M, r, du, n);
		incrementU(du, u, n);

		// spectral radius
		normDum1 = normDu;
		normDu = norm(du, n);
		printf("sp-rad: %4.2f \n", normDu / normDum1);

		// afficher_vect(u, n);
	}
	printf(" \n");
}
// project
void matrixVectorProduct(int n, const float A[MAX][MAX], const float in[MAX], float out[MAX])
{
	for (int i = 0; i < n; i++)
	{
		out[i] = 0;
		for (int j = 0; j < n; j++)
		{
			out[i] += A[i][j] * in[j];
		}
	}
}
void vectorAddition(int n, const float in[MAX], const float in2[MAX], float out[MAX])
{
	for (int i = 0; i < n; i++)
	{
		out[i] = in[i] + in2[i];
	}
}
void setZero(int n, float b[MAX])
{
	for (int i = 0; i < n; i++)
	{
		b[i] = 0;
	}
}
void setIC(int n, float b[MAX])
{
	for (int i = 0; i < n; i++)
	{
		if (i < n / 2)
		{
			b[i] = 0;
		}
		else
		{
			b[i] = 1;
		}
	}
}
void constructC(int n, float A[MAX][MAX], float dt, float l)
{
	// float l = 1;
	float dx = l / (n + 1);

	A[0][0] = 1 - 2 * dt / dx / dx;
	A[0][1] = 2*dt / dx / dx; //neu
	// A[0][1] = dt / dx / dx; //dir
	for (int i = 1; i < n-1; i++)
	{
		A[i][i] = 1 - 2 * dt / dx / dx;
		A[i][i + 1] = dt / dx / dx;
		A[i][i - 1] = dt / dx / dx;
	}
	A[n - 1][n - 1] = 1 - 2 * dt / dx / dx;
	A[n - 1][n - 2] = 2 * dt / dx / dx; //neu
	// A[n - 1][n - 2] = 1 * dt / dx / dx; //dir
}
void constructB(int n, float b[MAX], float dt, float dx)
{
	setZero(n, b);
	// b[0] = 0*dt/dx/dx; //dir
	b[0] = -2*2*dt/dx; //neu
	// b[n - 1] = 2*dt/dx/dx; //dir
	b[n - 1] = 2*2*dt/dx; //neu
}

void constructX(int n, float x[MAX], float xMin, float xMax)
{
	float dx = (xMax - xMin) / (n + 1);
	for (int i = 0; i < n; i++)
	{
		x[i] = xMin + dx + i * dx;
	}
}
int main()
{
	int n = 19;
	float dt = 1e-1;
	float xMin = -10;
	float xMax = 10;
	float dx = (xMax - xMin) / (n + 1);
	//	float x[MAX],b[MAX],y[MAX];
	float u[MAX], b[MAX], x[MAX];
	float C[MAX][MAX];
	// setZero(n,u);
	setIC(n, u);

	constructC(n, C, dt, xMax - xMin);
	constructB(n, b,dt,dx);
	constructX(n, x, xMin, xMax);

	int nt = 40/dt;

	for (int it = 0; it < nt; it++)
	{
		float work1[MAX], work2[MAX];
		matrixVectorProduct(n, C, u, work1);
		vectorAddition(n, work1, b, u);
		printf("t: %4.4f \t iter: %d \t D: %4.6f \t norm: %4.4f \n", it * dt+dt, it, dt / dx / dx, norm(u, n) * sqrtf(dx));
	}

	// save result
	FILE *fptr;
	fptr = fopen("data.dat", "w");
	for (int i = 0; i < n; i++)
	{
		fprintf(fptr, "%1.17e\t\t", x[i]);
		fprintf(fptr, "%1.17e\n", u[i]);
	}

	fclose(fptr);
}
int mainTP()
{
	int n;
	//	float x[MAX],b[MAX],y[MAX];
	float x[MAX], y[MAX];

	float L[MAX][MAX], U[MAX][MAX];
	//	float A[MAX][MAX],L[MAX][MAX],U[MAX][MAX];

	printf("Veuillez entrer la valeur de n : ");
	//	if (scanf("%d",&n) !=1) printf("error");
	n = 3;
	// float A[MAX][MAX] = {{1, 2, 3}, {2, 5, 10}, {3, 10, 26}};
	// float A[MAX][MAX] = {{0.7, -0.4}, {-0.2, 0.5}};
	//	float A[MAX][MAX]={{1,2,3},{0,1,4},{0,0,1}};
	float A[MAX][MAX] = {{4, -1, 0}, {-1, 4, -1}, {0, -2, 4}}; // example 3.1
	//	float A[MAX][MAX]={{1,0,0},{2,1,0},{3,4,1}};

	// float b[MAX] = {5, 6, 7};
	// float b[MAX] = {0.3, 0.3};
	float b[MAX] = {3.0 / 16, 4.0 / 16, 6.0 / 16}; // eg 3.1

	printf("Veuillez remplir la matrice A : \n");
	//	remplir_mat_2(A , n);
	//	remplir_mat_rand(A , n);

	afficher_mat(A, n);

	printf("Veuillez remplir le vecteur b : \n");
	//	remplir_vect_2(b,n);
	afficher_vect(b, n);

	factoriser_LU(A, L, U, n);
	printf("Matrice L : \n");
	afficher_mat(L, n);
	printf("Matrice U : \n");
	afficher_mat(U, n);

	printf("Solution de Ly=b: \n");
	//	resol_trig_inf( A,y,b ,n);
	resol_trig_inf(L, y, b, n);
	afficher_vect(y, n);

	printf("Solution de Ux=y: \n");
	//	resol_trig_sup( A,x,b ,n);
	resol_trig_sup(U, x, y, n);
	afficher_vect(x, n);

	//*/

	printf("resid \n");
	float r[MAX];
	calculateResid(A, b, x, r, n);
	// afficher_vect(r, n);
	printf("%4.2f \n", norm(r, n));

	// tp2
	float u2[MAX] = {0, 0, 0};
	solveIter(A, b, u2, n, 1e-4, 100);
	// afficher_vect(u2, n);
	return 0;
}
