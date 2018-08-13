#include <stdio.h>
#include <mpi.h>
#include <time.h>

#define IMAX 4
#define JMAX 4

//function to calculate output matrix
/*int output(int* A, int i, int j, int* h, int s){
  char temp = 0
  int p, q;
  for(p = -s; p < s; p++)
   for(q = -s; q < s; q++)
     temp += A[i-p,j-q];
  temp *= h[p+1,q+1];
  return temp;
}
*/
unsigned char** randMatr(int n,int m,int l){
  int i, j;
	time_t t;
	srand(time(NULL));

	unsigned char **A = malloc((n+2)*sizeof(unsigned char*));
	for(i = 0; i < n+2; i++)
		A[i] = malloc((m+2)*sizeof(char));

	for(i = 0; i < n+2; i++)
		for(j = 0; j < m+2; j++)
			if(i == 0 || j == 0 || i == n+1 || j == m+1)
				A[i][j] = l;
			else
				A[i][j] = l;//rand()%256;

	/*for(i = 0; i < n+2; i++){
		for(j = 0; j < m+2; j++)
			printf("%d,", A[i][j]);
		printf("\n");
	}*/
  return A;
}

void print_hallos(unsigned char** a,int r)
{
  printf("\n\n PROC %d\n\n",r );
  printf("first row:\n");
  for(int j=0;j<JMAX+2;j++) printf("%u ",a[0][j] );
  printf("\n\n first col\n\n");
  for(int i=0;i<IMAX+2;i++) printf("%u ",a[i][0] );
  printf("\n\n last row\n" );
  for(int j=0;j<JMAX+2;j++) printf("%u ",a[IMAX+1][j] );
  printf("\n\nlast col\n" );
  for(int i=0;i<IMAX+2;i++) printf("%u ",a[i][JMAX+1] );
  printf("\n\n END PROC %d\n",r );
}


int main(int argc,char** argv) {
   int my_rank, comm_sz;
   MPI_Datatype Row,Column;
   MPI_Status  status;

   //char a[IMAX][JMAX];


   MPI_Init(NULL, NULL);

   /* Derived Datatypes : Row & Column */
   MPI_Type_vector( IMAX, 1, JMAX+1, MPI_CHAR,&Column);     // +2 for the two halo columns and rows
   MPI_Type_contiguous(JMAX,MPI_CHAR,&Row);
   MPI_Type_commit(&Column);
   MPI_Type_commit(&Row);


   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   MPI_Comm new;
   int dims[2] = {4,4} ;
   int periods[2] = {0,0};
   MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,0,&new);

   unsigned char** a = randMatr(IMAX,JMAX,my_rank);

   int E,W,N,S,NE,NW,SE,SW;
   MPI_Cart_shift(new,0,1,&N,&S);
   MPI_Cart_shift(new,1,1,&W,&E);

   NE=N+1; NW=N-1; SE=S+1; SW=S-1;

   /* -1 same as MPI_Proc_NULL */
   if(N == -1)  NE=NW=-1;
   if(S == -1)  SE=SW=-1;
   if(E == -1)  NE=SE=-1;
   if(W == -1)  NW=SW=-1;


   MPI_Request reqsend[8], reqrecv[8];

   MPI_Isend(&a[1][1],1,Row,N,0,new,&reqsend[0]);       // Send row north
   MPI_Isend(&a[IMAX-1][1],1,Row,S,0,new,&reqsend[1]);    // Send row south
   MPI_Isend(&a[1][1],1,Column,W,0,new,&reqsend[2]);    // Send column west
   MPI_Isend(&a[1][JMAX-1],1,Column,E,0,new,&reqsend[3]); // Send column east
   MPI_Isend(&a[1][1],1,MPI_CHAR,NW,0,new,&reqsend[4]);
   MPI_Isend(&a[1][JMAX-1],1,MPI_CHAR,NE,0,new,&reqsend[5]);
   MPI_Isend(&a[IMAX-1][1],1,MPI_CHAR,SW,0,new,&reqsend[6]);
   MPI_Isend(&a[IMAX-1][JMAX-1],1,MPI_CHAR,SE,0,new,&reqsend[7]);


   MPI_Irecv(&a[0][1], 1, Row, N, 0, new, &reqrecv[0]);  // recv from north
   MPI_Irecv(&a[IMAX+1][1], 1, Row, S, 0, new, &reqrecv[1]); // recv from south
   MPI_Irecv(&a[1][0], 1, Column, W, 0, new, &reqrecv[2]); // recv from west
   MPI_Irecv(&a[1][JMAX+1], 1, Column, E, 0, new, &reqrecv[3]);  // recv from east
   MPI_Irecv(&a[0][0], 1, MPI_CHAR, NW, 0, new, &reqrecv[4]);
   MPI_Irecv(&a[0][JMAX+1], 1, MPI_CHAR, NE, 0, new, &reqrecv[5]);
   MPI_Irecv(&a[IMAX+1][0], 1, MPI_CHAR, SW, 0, new, &reqrecv[6]);
   MPI_Irecv(&a[IMAX+1][JMAX+1], 1, MPI_CHAR, SE, 0, new, &reqrecv[7]);


   // h' matrix
  /* float h[3][3] = {{0.0625, 0.125, 0.0625}, {0.125, 0.25, 0.125}, {0.0625, 0.125, 0.0625}};
   char** b = malloc((IMAX+1)*sizeof(char*));
   for(int k = 0; k < JMAX+1; k++)
      b[k] = malloc(JMAX+1);
   for(int i = 1; i <= IMAX-1; i++)
      for(int j = 1; j <= JMAX-1; j++)
          b[i,j] = output(a,i,j,h,1);*/



    MPI_Waitall(8,reqrecv,MPI_STATUSES_IGNORE);
    MPI_Waitall(8,reqsend,MPI_STATUSES_IGNORE);
    print_hallos(a,my_rank);


   MPI_Finalize();
   return 0;
}  /* main */
