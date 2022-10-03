#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 700
#include <time.h>





void factoriser_LU( const  float A[MAX][MAX] , float L[MAX][MAX] , float U[MAX][MAX], int n) 
{
//	int i,j,k;
	// factorisation à faire
for (int j=0; j<n; j++)
{
for (int i=0; i<j+1; i++)
{
U[i][j]=A[i][j];
//printf("%4.2f \t ",U[i][j]);
//printf("%4.2f \n ",A[i][j]);
for (int k=0; k<i; k++)
{
U[i][j]-=L[i][k]*U[k][j];
}//k
}//i
L[j][j]=1;
for (int i=j+1; i<n; i++)
{
L[i][j]=A[i][j];
//printf("%4.2f \n ",L[i][j]);

for (int k=0; k<j; k++)
{
L[i][j]-=L[i][k]*U[k][j];
}//k
L[i][j]=L[i][j]/U[j][j];
}//i

}//j
	
}


void resol_trig_inf( float A[MAX][MAX] , float x[MAX], float b[MAX],int n) {
//	int i,k;
	
	//la descente à faire
for (int i=0; i<n; i++)
{
x[i]=b[i];
for (int j=0; j<i; j++)
{
x[i]-=x[j]*A[i][j];
}
x[i]=x[i]/A[i][i];
}
}


void resol_trig_sup( float A[MAX][MAX] , float x[MAX], float b[MAX],int n) {
	// la remontée à faire
for (int i=n-1; i>-1; i--)
{
x[i]=b[i];
for (int j=n-1; j>i; j--)
{
x[i]-=x[j]*A[i][j];
}
x[i]=x[i]/A[i][i];
}

}	

void afficher_vect( float x[MAX],int n) {
	int i;
	for (i=0;i<n;i++) printf("%4.2f\n",x[i]);
}


void afficher_mat(float A[MAX][MAX] , int n) {
	int i,j;
	
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
			printf("%.1f\t",A[i][j]);
		}
		printf("\n");
	}
}

void remplir_vect( float x[MAX],int n) {
	int i;
	for (i=0;i<n;i++) {
		printf("valeur [%d] : ",i+1);
	    if (scanf("%f",&x[i]) !=1) printf("error");
		
	}
		
}

void remplir_mat(float A[MAX][MAX] , int n) {
	int i,j;

	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
		printf("valeur [%d][%d] : ",i+1,j+1);
	    if (scanf("%f",&A[i][j]) !=1) printf("error");
		}
		printf("\n");
	}
}

void remplir_mat_rand(float A[MAX][MAX] , int n) {

srand(time(NULL));
	int i,j;

	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
float r=rand();
A[i][j]=r/RAND_MAX;
		}
	}
}




int main() {
	int n;
//	float x[MAX],b[MAX],y[MAX];
	float x[MAX],y[MAX];

	float L[MAX][MAX],U[MAX][MAX];
//	float A[MAX][MAX],L[MAX][MAX],U[MAX][MAX];
	
	printf("Veuillez entrer la valeur de n : ");
//	if (scanf("%d",&n) !=1) printf("error");
	n=3;
	float A[MAX][MAX]={{1,2,3},{2,5,10},{3,10,26}};
//	float A[MAX][MAX]={{1,2,3},{0,1,4},{0,0,1}};
//	float A[MAX][MAX]={{1,0,0},{2,1,0},{3,4,1}};

	float b[MAX]={5,6,7};

	printf("Veuillez remplir la matrice A : \n");
//	remplir_mat_2(A , n);
//	remplir_mat_rand(A , n);

	afficher_mat(A,n);
		
	printf("Veuillez remplir le vecteur b : \n");
//	remplir_vect_2(b,n);
	afficher_vect(b,n);	
	

	factoriser_LU(A,L,U,n);
	printf("Matrice L : \n");
	afficher_mat(L , n);
	printf("Matrice U : \n");
	afficher_mat(U , n);
	
	printf("Solution de Ly=b: \n");
//	resol_trig_inf( A,y,b ,n);
	resol_trig_inf( L,y,b ,n);
	afficher_vect(y,n);	
	
	
	printf("Solution de Ux=y: \n");
//	resol_trig_sup( A,x,b ,n);
	resol_trig_sup( U,x,y ,n);
	afficher_vect(x,n);
	
	//*/
	



	return 0;
}
