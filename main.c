#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M 3 //lines
#define N 4 //

void printMatrix(int m, int n,float* matrix) 
{
	int i, j;
    for(i = 0; i < m; i++)
	{
        for(j = 0; j < n; j++)
		{
            printf(" %f\t", *(matrix + (m+1)*i + j));
        }
        printf("\n");
    } 
    printf("\n");
}

void modMatrixLine(int m, int n, float* matrix, float modifier, int indx) // 2 - 1 
{
	int i, j;
	for(j = 0; j < n; j++)
	{
		*(matrix + (m+1)*indx + j) = modifier*(*(matrix + (m+1)*indx + j));
	}
}

void SumMatrixLine(int m, int n, float* matrix, float k1, float k2, int indx1, int indx2) // 2 - 1 
{
	int i, j;
	for(j = 0; j < n; j++)
	{
		*(matrix + (m+1)*indx2 + j) = k2*(*(matrix + (m+1)*indx2 + j)) + k1*(*(matrix + (m+1)*indx1 + j));
	}
}

int main(void) 
{
    int i, j, maxIndex;
	float matrix[M][N] = {0};
    float max, buff;
    matrix[0][0] = 0.14; 	matrix[0][1] = 0.24; 	matrix[0][2] = -0.84; 	matrix[0][3] = 1.11; //0.14*x + 0.24*y - 0.84*z = 1.11
	matrix[1][0] = 1.07; 	matrix[1][1] = -0.83; 	matrix[1][2] = 0.56; 	matrix[1][3] = 0.48; //1.07*x - 0.83 + 0.56*z = -0.84
	matrix[2][0] = 0.64; 	matrix[2][1] = 0.43;	matrix[2][2] = -0.38;	matrix[2][3] = -0.83; //0.64*x + 0.43*y - 0.38*z = 0.64
    printMatrix (M,N,(float*) matrix);

    for(i = 0; i < M; i++)
	{
		for(j = 0; j < N; j++)
        {
            if(j = 0)
            {
                max = matrix[i][0];
                continue;
            }
            if(max < matrix[i][j])
            {
                max = matrix[i][j];
                maxIndex = j;
            }
        }
        if(maxIndex = i)
        {
            continue;
        }
        buff = matrix[i][i];
        matrix[i][i] = matrix[i][maxIndex];
        matrix[i][maxIndex] = buff;
	}

    printMatrix (M,N,(float*) matrix);

    for(i = 0; i < M; i++)
	{
		modMatrixLine(M, N, (float*) matrix, 1.0/matrix[i][i], i);
	}

    printMatrix (M,N,(float*) matrix);

    /*
	int i, j;
	float matrix[M][N] = {0};
	
//	matrix[0][0] = 1; 	matrix[0][1] = 1; 	matrix[0][2] = 2;
//	matrix[1][0] = 2; 	matrix[1][1] = 4; 	matrix[1][2] = 6;

	matrix[0][0] = 0.14; 	matrix[0][1] = 0.24; 	matrix[0][2] = -0.84; 	matrix[0][3] = 1.11; //0.14*x + 0.24*y - 0.84*z = 1.11
	matrix[1][0] = 1.07; 	matrix[1][1] = -0.83; 	matrix[1][2] = 0.56; 	matrix[1][3] = 0.48; //1.07*x - 0.83 + 0.56*z = -0.84
	matrix[2][0] = 0.64; 	matrix[2][1] = 0.43;	matrix[2][2] = -0.38;	matrix[2][3] = -0.83; //0.64*x + 0.43*y - 0.38*z = 0.64
	
	
	printMatrix (M,N,(float*) matrix);
	
	for(j = 0; j < N-1; j++)
	{
		modMatrixLine(M, N, (float*) matrix, 1.0/matrix[j][j], j);
		printMatrix (M,N,(float*) matrix);
		for(i = 0; i < M - 1 - j; i++)
		{
			SumMatrixLine(M, N, (float*) matrix, - matrix[i+1+j][j], 1.0, j, i+1+j);
			printMatrix (M,N,(float*) matrix);
		}
	}
	
	float x[M] = {0};
	
	for(i = (M - 1); i >= 0; i--)
	{
		x[i] = matrix[i][N-1];
	//	printf("%f\n", x[i]);
		for(j = 0; j < N - 2 - i; j++)//j = i; j < N - 2; j++)
		{
			if(i == M-1)
			{
				continue;
			}
		//	printf("%f,  %f\n", x[i+j+1],matrix[i][i+j+1]);
			x[i] = x[i] - x[i+j+1]*matrix[i][i+j+1];
		}
	}	
	
	printf("\n");
	printf("\n");
	printMatrix(1,M,x);

	
	*/
	return 0;
}
