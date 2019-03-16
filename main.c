#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M 3 //lines
#define N 4 //

#define ACCURASY 0.00001 //

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

void MovMatrixLine(int m, int n, float* matrix, int indx1, int indx2) // 2 - 1 
{
	int i, j;
	float buff;
	for(j = 0; j < n; j++)
	{
		buff = *(matrix + (m+1)*indx2 + j);
		*(matrix + (m+1)*indx2 + j) = (*(matrix + (m+1)*indx1 + j));
		(*(matrix + (m+1)*indx1 + j)) = buff;
	}
}

void SimpleIterationMode()
{
	 int i, j, maxIndex;
	float matrix[M][N] = {};
	float x[M] = {};
	float oldx[M] = {};
    float max, buff;
    
    matrix[0][0] = 0.12; 	matrix[0][1] = -0.43; 	matrix[0][2] = 0.14; 	matrix[0][3] = -0.17; //0.14*x + 0.24*y - 0.84*z = 1.11
	matrix[1][0] = -0.07; 	matrix[1][1] = 0.34; 	matrix[1][2] = 0.72; 	matrix[1][3] = 0.62; //1.07*x - 0.83 + 0.56*z = -0.84
	matrix[2][0] = 1.18; 	matrix[2][1] = -0.08;	matrix[2][2] = -0.25;	matrix[2][3] = 1.12; //0.64*x + 0.43*y - 0.38*z = 0.64

    for(i = 0; i < M; i++)
	{
		for(j = 0; j < N - 1; j++)
        {
            if(j == 0)
            {
                max = matrix[i][0];
            	maxIndex = 0;
                continue;
            }
            if(fabs(max) < fabs(matrix[i][j]))
            {
                max = matrix[i][j];
                maxIndex = j;
            }
        }
        if(i == maxIndex)
        {
        	continue;
		}
		else
		{
			MovMatrixLine(M,N,(float*) matrix,i,maxIndex);
		}
	}

    for(i = 0; i < M; i++)
	{
		modMatrixLine(M, N, (float*) matrix, -1.0/matrix[i][i], i);
		matrix[i][i] = 0;
	}
	
	
	for(i = 0; i < M; i++)
	{
		x[i] = 0;
		oldx[i] = 0;
	}
	
    float delta;
    int k = 0;
    do
    {
    	for(i = 0; i < M; i++)
    	{
    		oldx[i] = x[i];
		}
    	for(i = 0; i < M; i++)
    	{
    		x[i] = 0;
    		for(j = 0; j < M; j++)
    		{
    			x[i] += matrix[i][j] * oldx[j];
			}
			x[i] -= matrix[i][N-1];
		}
		
		
		delta = fabs(x[0] - oldx[0]);
		for(i = 1; i < M; i++)
		{
			if(delta < fabs(x[i] - oldx[i]))
			{
				delta = fabs(x[i] - oldx[i]);
			}
		}
		k++;
	}
	while((fabs(delta) > ACCURASY));
//	printf("iteration counter = %d\n",k);
	printMatrix(1,M,x);
}

void ZeydelIterationMode()
{
	 int i, j, maxIndex;
	float matrix[M][N] = {};
	float x[M] = {};
	float oldx[M] = {};
    float max, buff;
    
    matrix[0][0] = 0.12; 	matrix[0][1] = -0.43; 	matrix[0][2] = 0.14; 	matrix[0][3] = -0.17; //0.14*x + 0.24*y - 0.84*z = 1.11
	matrix[1][0] = -0.07; 	matrix[1][1] = 0.34; 	matrix[1][2] = 0.72; 	matrix[1][3] = 0.62; //1.07*x - 0.83 + 0.56*z = -0.84
	matrix[2][0] = 1.18; 	matrix[2][1] = -0.08;	matrix[2][2] = -0.25;	matrix[2][3] = 1.12; //0.64*x + 0.43*y - 0.38*z = 0.64

    for(i = 0; i < M; i++)
	{
		for(j = 0; j < N - 1; j++)
        {
            if(j == 0)
            {
                max = matrix[i][0];
            	maxIndex = 0;
                continue;
            }
            if(fabs(max) < fabs(matrix[i][j]))
            {
                max = matrix[i][j];
                maxIndex = j;
            }
        }
        if(i == maxIndex)
        {
        	continue;
		}
		else
		{
			MovMatrixLine(M,N,(float*) matrix,i,maxIndex);
		}
	}

    for(i = 0; i < M; i++)
	{
		modMatrixLine(M, N, (float*) matrix, -1.0/matrix[i][i], i);
		matrix[i][i] = 0;
	}
	
	
	for(i = 0; i < M; i++)
	{
		x[i] = 0;
		oldx[i] = 0;
	}
	
    float delta;
    int k = 0;
    do
    {
    	for(i = 0; i < M; i++)
    	{
    		oldx[i] = x[i];	
		}
    	for(i = 0; i < M; i++)
    	{
    		x[i] = 0;
    		for(j = 0; j < M; j++)
    		{
    			x[i] += matrix[i][j] * oldx[j];
			}
			x[i] -= matrix[i][N-1];
			oldx[i] = x[i];
		}
		
		
		delta = fabs(x[0] - oldx[0]);
		for(i = 1; i < M; i++)
		{
			if(delta < fabs(x[i] - oldx[i]))
			{
				delta = fabs(x[i] - oldx[i]);
			}
		}
		k++;
	}
	while((fabs(delta) > ACCURASY));
//	printf("iteration counter = %d\n",k);
	printMatrix(1,M,x);
}

int main(void) 
{
   	printf("Reshit sistemy:\n1 - metodom postyh iteraciy\n2 - metodom Zeydelya\n");
   	char choice = getch();
 	if(choice = '1')
 	{
 		SimpleIterationMode();
	}
	else
	{
		ZeydelIterationMode();
	}
	return 0;
}
