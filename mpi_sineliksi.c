#include <stdio.h>
#include <mpi.h>
#include <time.h>

#define IMAX 4
#define JMAX 4

//function to calculate output matrix
int output(int* A, int i, int j, int* h, int s){
  char temp = 0
  int p, q;
  for(p = -s; p < s; p++)
   for(q = -s; q < s; q++)
     temp += A[i-p,j-q];
  temp *= h[p+1,q+1];
  return temp;
}

char** randMatr(int n,int m){
  int i, j;
	time_t t;
	srand(time(NULL));

	char **A = malloc((n+2)*sizeof(char*));
	for(i = 0; i < n+2; i++)
		A[i] = malloc((m+2)*sizeof(char));

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
  return A;
}

int main(int argc,char** argv) {
   int my_rank, comm_sz;
   MPI_Datatype Row,Column;
   MPI_Status  status;

   char a[IMAX][JMAX];
   char** test = randMatr(IMAX,JMAX);

   MPI_Init(NULL, NULL);

   /* Derived Datatypes : Row & Column */
   MPI_Type_vector( IMAX, 1, JMAX, MPI_CHAR,&Column);
   MPI_Type_contiguous(JMAX,MPI_CHAR,&Row);
   MPI_Type_commit(&Column);
   MPI_Type_commit(&Row);


   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


   int E,W,N,S,NE,NW,SE,SW;
   E=my_rank+1; W=my_rank-1; N=my_rank-JMAX; S=my_rank+JMAX;
   NE=N+1; NW=N-1; SE=S+1; SW=S-1;

   if(my_rank<=JMAX){
     N=NE=NW=MPI_PROC_NULL;
   }
   if(my_rank>=comm_sz-JMAX){
     S=SE=SW=MPI_PROC_NULL;
   }
   if(my_rank%JMAX == 0){
     W=SW=NW=MPI_PROC_NULL;
   }
   if(my_rank%JMAX == JMAX-1){
      E=SE=NE=MPI_PROC_NULL;
   }

   MPI_Request req;

   MPI_Isend(&a[1][1],1,Row,N,0,MPI_COMM_WORLD,&req);       // Send row north
   MPI_Isend(&a[IMAX-1][1],1,Row,S,0,MPI_COMM_WORLD,&req);    // Send row south
   MPI_Isend(&a[1][1],1,Column,W,0,MPI_COMM_WORLD,&req);    // Send column west
   MPI_Isend(&a[1][JMAX-1],1,Column,E,0,MPI_COMM_WORLD,&req); // Send column east
   MPI_Isend(&a[1][1],1,MPI_CHAR,NW,0,MPI_COMM_WORLD,&req);
   MPI_Isend(&a[1][JMAX-1],1,MPI_CHAR,NE,0,MPI_COMM_WORLD,&req);
   MPI_Isend(&a[IMAX-1][1],1,MPI_CHAR,SW,0,MPI_COMM_WORLD,&req);
   MPI_Isend(&a[IMAX-1][JMAX-1],1,MPI_CHAR,SE,0,MPI_COMM_WORLD,&req);

   MPI_Irecv(&a[0][1], 1, Row, N, 0, MPI_COMM_WORLD, &req);  // recv from north
   MPI_Irecv(&a[IMAX][1], 1, Row, S, 0, MPI_COMM_WORLD, &req); // recv from south
   MPI_Irecv(&a[1][0], 1, Column, W, 0, MPI_COMM_WORLD, &req); // recv from west
   MPI_Irecv(&a[1][JMAX], 1, Column, E, 0, MPI_COMM_WORLD, &req);  // recv from east
   MPI_Irecv(&a[0][0], 1, MPI_CHAR, NW, 0, MPI_COMM_WORLD, &req);
   MPI_Irecv(&a[0][JMAX], 1, MPI_CHAR, NE, 0, MPI_COMM_WORLD, &req);
   MPI_Irecv(&a[IMAX][0], 1, MPI_CHAR, SW, 0, MPI_COMM_WORLD, &req);
   MPI_Irecv(&a[IMAX][JMAX], 1, MPI_CHAR, SE, 0, MPI_COMM_WORLD, &req);

   // h' matrix
   float h[3][3] = {{0.0625, 0.125, 0.0625}, {0.125, 0.25, 0.125}, {0.0625, 0.125, 0.0625}};
   char** b = malloc((IMAX+1)*sizeof(char*));
   for(int k = 0; k < JMAX+1; k++)
      b[k] = malloc(JMAX+1);
   for(int i = 1; i <= IMAX-1; i++)
      for(int j = 1; j <= JMAX-1; j++)
          b[i,j] = output(a,i,j,h,1);


   MPI_Finalize();
   return 0;
}  /* main */
