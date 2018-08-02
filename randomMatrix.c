#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void foo(int n, int m)
{
	int i, j;
	time_t t;
	srand(time(&t));
	
	int **A = malloc((n+2)*sizeof(int*));
	for(i = 0; i < n+2; i++)
		A[i] = malloc((m+2)*sizeof(int));
		
	for(i = 0; i < n+2; i++)
		for(j = 0; j < m+2; j++)
			if(i == 0 || j == 0 || i == n+1 || j == m+1)
				A[i][j] = 1;
			else
				A[i][j] = rand()%256;
				
	for(i = 0; i < n+2; i++){
		for(j = 0; j < m+2; j++)
			printf("%d,", A[i][j]);
		printf("\n");
	}
	
}

int main(void)
{
	foo(5,5);
}
