/* File:     mpi_output.c
 *
 * Purpose:  A program in which multiple MPI processes try to print
 *           a message.
 *
 * Compile:  mpicc -g -Wall -o mpi_output mpi_output.c
 * Usage:    mpiexec -n<number of processes> ./mpi_output
 *
 * Input:    None
 * Output:   A message from each process
 *
 * IPP:      Section 3.3.1  (pp. 97 and ff.)
 */
#include <stdio.h>
#include <mpi.h>
#include <time.h>

#define IMAX 25
#define JMAX 12

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


// SEND_PROTOTYPE: src=buffer_to_send, n=elements_number_to_send, rank=rank_of_dest
// RECV_PROTOTYPE: dest=buffer_to_recv, n=elements_number_to_recv, rank=rank_of_src

//=====================================================================================
void send_east(char** src,int dest)
{
  MPI_Request rq;
  MPI_Isend(&src[1][JMAX-1],1,Column,dest,0,MPI_COMM_WORLD,&rq);
}

void recv_east(char** arr,int source)
{
  MPI_Irecv(&arr[1][JMAX-1],1,Column,0,MPI_COMM_WORLD,&rq);
}

void send_north(int* src, int n, int rank){
  MPI_Isend(s, n, MPI_INT, rank, 0, MPI_COMM_WORLD);
}

void send_south(int* src, int n, int rank){
  MPI_Send(src, n, MPI_INT, rank, 0, MPI_COMM_WORLD);
}
void recv_south(int* dest, int n, int rank, MPI_Status status){
  MPI_Recv(dest, n, MPI_INT, rank, 0, MPI_COMM_WORLD , &status);
  MPI_Barrier(MPI_COMM_WORLD);
}

//=====================================================================================

void send_southwest(int* src, int rank){
  MPI_Send(src, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
}
void recv_southwest(int* dest, int rank, MPI_Status status){
  MPI_Recv(dest, 1, MPI_INT, rank, 0, MPI_COMM_WORLD , &status);
  MPI_Barrier(MPI_COMM_WORLD);
}

//====================================================================================

void send_southeast(int* src, int rank){
  MPI_Send(src, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
}
void recv_southeast(int* dest, int rank, MPI_Status status){
  MPI_Recv(dest, 1, MPI_INT, rank, 0, MPI_COMM_WORLD , &status);
  MPI_Barrier(MPI_COMM_WORLD);
}

//=====================================================================================

void send_west(int* src, int n, int rank){
  MPI_Send(src, n, MPI_INT, rank, 0, MPI_COMM_WORLD);
}
void recv_west(int* dest, int n, int rank, MPI_Status status){
  MPI_Recv(dest, n, MPI_INT, rank, 0, MPI_COMM_WORLD , &status);
  MPI_Barrier(MPI_COMM_WORLD);
}

//=====================================================================================

int main(void) {
   int my_rank, comm_sz;
   MPI_Datatype Row,Column;
   MPI_Status  status;

   char a[IMAX][JMAX];
   char** test = randMatr(IMAX,JMAX);

   for(int i=0;i<IMAX;i++){
     for(int j=0;j<JMAX;j++){
       a[i][j] = 9;
     }

   }

   MPI_Init(NULL, NULL);

   MPI_Type_vector( IMAX, 1, JMAX, MPI_CHAR,&Column);
   MPI_Type_contiguous(JMAX,MPI_CHAR,&Row);
   MPI_Type_commit(&Column);
   MPI_Type_commit(&Row);




   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   printf("Proc %d of %d > Does anyone have a toothpick?\n",
         my_rank, comm_sz);

   if(my_rank != 0){
     int k=1997;
     MPI_Send(&k, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
     MPI_Send(a, 1, Column, 0, 0, MPI_COMM_WORLD);
     printf("send ok\n" );
   }  else {

     char dat[IMAX][JMAX];
     int l=0;
     printf("***********************\n" );
     MPI_Recv(&l, 1, MPI_INT, 1, 0,MPI_COMM_WORLD , &status);
     printf("-------> L = %d\n",l );
     MPI_Recv(dat, 1, Column, 1, 0,MPI_COMM_WORLD , &status);
     printf("+++++++++++++++++++++++++\n"
      );
     for(int i=0;i<JMAX;i++){
       printf("dat[%d]=%d\n",i,dat[i][0] );
     }
   }


   int E,W,N,S,NE,NW,SE,SW;
   E=my_rank+1; W=my_rank-1; N=my_rank-JMAX; S=my_rank+JMAX;
   NE=N+1; NW=N-1; SE=S+1; SW=S-1;

   // to SRC tha einai o pointer sto antikeimeno pou theloume. Tha to dinoume san orisma apeytheias
   // to N tha einai o arithmos twn antikeimenwn analoga me thn periptwsh

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


   MPI_Finalize();
   return 0;
}  /* main */
