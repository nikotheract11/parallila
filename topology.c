#include <stdio.h>
#include <mpi.h>
#include <time.h>



int main(int argc,char** argv) {
   int my_rank, comm_sz;
   MPI_Datatype Row,Column;
   MPI_Status  status;

   //char a[IMAX][JMAX];


   MPI_Init(NULL, NULL);


   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   MPI_Comm new;
   int dims[2] = {4,4} ;
   int periods[2] = {0,0};
   MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,0,&new);

   int coords[2] = {0,2};
   int rank;
   if(my_rank==1){
     int s,d;
     MPI_Cart_shift(new,0,2,&s,&d);
     printf("proc:%d source = %d and dest=%d\n",my_rank,s,d );
   }

   //int a[3][3]={{1,2,3},{1,2,3},{1,2,3}};
   int **a = malloc(3*sizeof(int*));
   for(int i=0;i<3;i++) {
     a[i] = malloc(3*sizeof(int));
     for(int j=0;j<3;j++) a[i][j] = j;
   }

   MPI_Type_vector( 3, 1, 3, MPI_INT,&Column);     // +2 for the two halo columns and rows
   MPI_Type_commit(&Column);
   int **b = malloc(3*sizeof(int*));
   for(int i=0;i<3;i++) {
     b[i] = malloc(3*sizeof(int));
  //   for(int j=0;j<3;j++) a[i][j] = j;
   }
   if(my_rank==1) {
     MPI_Send(&a[0][2],1,Column,2,0,new);
   }
   if(my_rank==2){
     MPI_Recv(b[0]+1,3,Column,1,0,new,MPI_STATUSES_IGNORE);
     for(int o=0;o<3;o++) printf("%d\n",b[o][1] );
   }


      MPI_Finalize();
      return 0;
   }
