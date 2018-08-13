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


      MPI_Finalize();
      return 0;
   }
