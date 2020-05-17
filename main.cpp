#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "winbgi2.h"
#include "gauss.h"

void HilbertMatrix(int N, double** H);
void computeMatrix(int N, double** K);
void displayMatrix(int N, double** H);
void computeVec(int N, double** H, double* b);
void computeVector(int N, double* F);
void plotVec(int N, double* v);
void scan(int* x);
void wykres(double *T, int N);
double wyznacznik(int N, double** K);

double Tp = 273;
double Tk = 300;
double lambda = 58.0;
double w1 = 1.0, w2 = -2.0, w3 = 1.0;			//wspolczynniki w rownaniu
double L = 1.0;

void main()
{
	int n;
	double** H, ** K;
	double* b, * x, * F, * T;

	printf("Zadanie 1. macierz Hilberta \n\n");
	printf("Podaj ilosc rownan [N]: \n");
	scan(&n);

	H = (double**)malloc(n * sizeof(double*));
	K = (double**)malloc((n+1) * sizeof(double*));

	b = (double*)malloc(n * sizeof(double));
	x = (double*)malloc(n * sizeof(double));
	F = (double*)malloc((n+1) * sizeof(double));
	T = (double*)malloc((n+1) * sizeof(double));

	if (H == NULL && b == NULL && x == NULL && F == NULL && T == NULL)
	{
		printf("Wystapil blad alokacji :( \n");
		exit(1);
	}
	
	for (int i = 0; i < n+1; i++)
	{
		if (i < n)
		{
			H[i] = (double*)malloc(n * sizeof(double));
		}
		K[i] = (double*)malloc((n+1) * sizeof(double));

		if (H[i] == NULL && K[i] == NULL)
		{
			printf("Wystapil blad alokacji :( \n");
			exit(1);
		}
	}


	HilbertMatrix(n, H);
	printf("MACIERZ HILBERTA: \n");
	displayMatrix(n, H);
	computeVec(n, H, b);
	printf("WEKTOR B: \n");
	plotVec(n, b);
	gauss(n, H, x, b);
	printf("WEKTOR X: \n");
	plotVec(n, x);


	printf("Zadanie 2. rozklad temperatury w precie \n\n");

	computeMatrix(n + 1, K);
	printf("MACIERZ K: \n");
	displayMatrix(n+1, K);
	computeVector(n, F);
	printf("WEKTOR F: \n");
	plotVec(n+1,F);
	gauss(n+1, K, T, F);
	printf("WEKTOR T: \n");
	plotVec(n+1, T);
	wykres(T, n);
	
	printf("Wyznacznik macierzy K[%d][%d] = %lf \n",n+1,n+1, wyznacznik(n+1, K));
	

	free(T);
	free(F);
	free(x);
	free(b);
	for (int k = 0; k < n+1; k++)
	{
		if (k < n)
		{
			free(H[k]);
		}
		free(K[k]);
	}
	free(H);
	free(K);

	wait();
}

void scan(int* x)		//funkcja do wczytywania i sprawdzania zgodnosci inputu
{

	int temp;
	while (scanf("%d", &temp) != 1)
	{
		printf("Nieprawid³owy format danych! \n");
		int n;
		while ((n = getchar()) != EOF && n != '\n');
	}
	*x = temp;

}

void HilbertMatrix(int N, double** H)
{

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			H[i][j] = 1./(1 + i + j);
		}
	}

}

void computeMatrix(int N, double** K)
{

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			K[i][j] = 0.0;
		}
	}
	for (int k=0; k<N-2;k++)
	{
		K[k + 1][k] =w1;
		K[k + 1][k + 1] = w2;
		K[k + 1][k + 2] = w3;
	}
	K[0][0] = 1.0;
	K[N - 1][N - 1] = 1.0;

}

void displayMatrix(int N, double** H)
{

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("[%lf] ", H[i][j]);
		}
		printf("\n");
	}
	printf("\n");

}

void computeVec(int N, double** H, double* b)
{

	for (int i = 0; i < N; i++)
	{
		b[i] = 0;
		for (int j = 0; j < N; j++)
		{
			b[i] += H[i][j];
		}
	}

}

void computeVector(int N, double* F)
{

	double h = L / N;
	double x = h;
	for (int i = 1; i < N; i++)
	{
		F[i] = ((-pow(10.0, 4.0) * sin(x * 3.1416)) / lambda )* pow(h, 2.0);
		x += h;
	}
	F[0] = Tp;
	F[N] = Tk;

}

void plotVec(int N, double* v)
{

	for (int i = 0; i < N; i++)
	{
		printf("[%lf]\n",v[i]);
	}
	printf("\n");

}

void wykres(double* T, int N)
{
	double x = 0;

	graphics(600, 400);
	scale(0, 27*10, 1, 31*10);

	for (int d = 0; d < N + 1; d++)
	{
		point(x, T[d]);
		circle(x, T[d], 5);
		lineto(x, T[d]);
		x += L / N;
	}

}

double wyznacznik(int N, double** K)
{

	double r, wyz = 1;
	for (int i = 0; i < N; i++) 
	{
		for (int j = 0; j < N; j++)
		{
			if (j > i) 
			{
				r = K[j][i] / K[i][i];
				for (int k = 0; k < N; k++) 
				{
					K[j][k] -= r * K[i][k];
				}
			}
		}
	}
	for (int i = 0; i < N; i++)
	{
		wyz *= K[i][i];
	}
	return wyz;

}