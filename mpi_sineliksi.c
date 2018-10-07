#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

//function to calculate output matrix
static inline int output(unsigned char* A, int i, int j, unsigned char** h, int s,int ROWS){
  unsigned char temp = 0;
  int p, q;
  for(p = -s; p < s; p++)
   for(q = -s; q < s; q++)
     temp += A[(i-p)*ROWS+j-q]*h[p+1][q+1];
  return temp;
}

void resize(char *file,int new_i,int new_j,int i,int j)
{
   FILE *fp;
   FILE *w;

   int off=0;
   int ofw=0;

   fp=fopen(file,"r");
   w=fopen("new.raw","w");

   char * buf = malloc(new_i * sizeof(char));

   for(int k=0;k<new_j;k++){

   	fread(buf,1,new_i,fp);
   	fwrite(buf,1,new_i,w);

      off += i;
      ofw += new_i;

  	   fseek(fp,off,SEEK_SET);
      fseek(w,ofw,SEEK_SET);
   }

   fclose(fp);
   fclose(w);

}

static inline char equals(unsigned char* before,unsigned char* after,int IMAX,int JMAX){
  for(int i=1;i<IMAX+1;i++){
    for(int j=1;j<JMAX+1;j++)    {
      if(after[i*(JMAX+2)+j] != before[i*(JMAX+2)+j]) return 0;
    }
  }
  return 1;
}

static inline void mat(unsigned char** a,int IMAX,int JMAX){
  time_t t;
	srand(time(NULL));

  *a = malloc((IMAX+2)*(JMAX+2)*sizeof(unsigned char));
  for(int i=0;i<IMAX+2;i++) {
    for(int j=0;j<JMAX+2;j++) {
      (*a)[i*(JMAX+1)+j]=rand()%256;
    }
  }
}

static inline unsigned char* fit(unsigned char *buf,int IMAX,int JMAX,int pix)
{
   unsigned char* tmp = malloc((IMAX+2)*(JMAX+2*pix)*sizeof(unsigned char));
   for(int i=0;i<IMAX+2;i++) {
     for(int j=0;j<JMAX+2*pix;j++) {
       if(i==0 || i==IMAX-1 || j==0 || j==JMAX-1)
        tmp[i*(JMAX+2*pix)+j] = 1;
     }
  }
   for(int i=1;i<IMAX+1;i++) {
     for(int j=1;j<JMAX+1*pix;j++) {
        tmp[i*(JMAX+2*pix)+j] = buf[(i-1)*(JMAX)+j];
     }
  }

  return tmp;
}

static inline unsigned char* repeat(unsigned char *buf,int IMAX,int JMAX,int pix)
{
   unsigned char* tmp = malloc(IMAX*JMAX*sizeof(unsigned char));
   for(int i=1;i<IMAX+1;i++) {
     for(int j=1;j<JMAX+1*pix;j++) {
        tmp[(i-1)*JMAX+j-1] = buf[i*(JMAX+2*pix)+j];
     }
  }

  return tmp;
}

int main(int argc,char** argv) {
   int my_rank, comm_sz;
   MPI_Datatype Row,Column;
   MPI_Status status;

   char* filenamein = argv[3];
   char* filenameout = "output.raw";

   resize(filenamein,1000,1200,1920,2520);

   int IMAX = atoi(argv[1]);
   int JMAX = atoi(argv[2]);
   int pix = atoi(argv[4]);
   int ROWS = JMAX+2;

   // calculate size of the file
   FILE *fp;
   fp=fopen(filenamein,"r");
   fseek(fp, 0L, SEEK_END);
   int FILESIZE = (int)ftell(fp);
   fclose(fp);


   MPI_Init(NULL, NULL);


   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   IMAX /= sqrt(comm_sz);
//   IMAX *= pix;
   JMAX /= sqrt(comm_sz);
   JMAX *= pix;
   ROWS = JMAX + 2*pix;

   MPI_Comm new;
   int dims[2] = {sqrt(comm_sz),sqrt(comm_sz)} ;
   int periods[2] = {0,0};
   MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,0,&new);

   /* Derived Datatypes : Row & Column */
   MPI_Type_vector(IMAX,1,JMAX+2*pix,MPI_UNSIGNED_CHAR,&Column);     // +2 for the two halo columns and rows
   MPI_Type_commit(&Column);
   MPI_Type_contiguous(JMAX,MPI_UNSIGNED_CHAR,&Row);
   MPI_Type_commit(&Row);

   unsigned char* a ;

   // open and read input image
   MPI_File f;
   MPI_File_open(new,filenamein,MPI_MODE_RDONLY,MPI_INFO_NULL,&f);
   int BUFSIZE = IMAX*JMAX;
   MPI_Comm_rank(new, &my_rank);

   int offsetb = 0;
   MPI_Offset offsetf = (my_rank/sqrt(comm_sz))*sqrt(comm_sz)*IMAX*JMAX + (my_rank%(int)sqrt(comm_sz))*JMAX;
   MPI_File_seek(f,offsetf,MPI_SEEK_SET);

   unsigned char* buffer = malloc(BUFSIZE*sizeof(unsigned char));

   for(int i = 0; i < IMAX; i++){
     MPI_File_read(f,&buffer[offsetb],JMAX,MPI_UNSIGNED_CHAR,&status);
     offsetf += JMAX*sqrt(comm_sz);
     offsetb += JMAX;
     MPI_File_seek(f,offsetf,MPI_SEEK_SET);
   }
   MPI_File_close(&f);

   /* a contains buffer and also halo points */
   a = fit(buffer,IMAX,JMAX,pix);

   int E,W,N,S,NE,NW,SE,SW;
   MPI_Cart_shift(new,0,1,&N,&S);
   MPI_Cart_shift(new,1,1,&W,&E);

   NE=N+1; NW=N-1; SE=S+1; SW=S-1;

   /* -1 same as MPI_Proc_NULL */
   if(N == MPI_PROC_NULL)  NE=NW=MPI_PROC_NULL;
   if(S == MPI_PROC_NULL)  SE=SW=MPI_PROC_NULL;
   if(E == MPI_PROC_NULL)  NE=SE=MPI_PROC_NULL;
   if(W == MPI_PROC_NULL)  NW=SW=MPI_PROC_NULL;


   MPI_Request reqsend[8], reqrecv[8];


  MPI_Barrier(new);
  double start = MPI_Wtime();
  char* b = malloc((IMAX+2*pix)*(JMAX+2*pix)*sizeof(unsigned char));


   while(1){

	   MPI_Isend(&a[1*ROWS+1*pix],1,Row,N,0,new,&reqsend[0]);       // Send row north
	   MPI_Isend(&a[IMAX*ROWS+1*pix],1,Row,S,0,new,&reqsend[1]);    // Send row south
	   MPI_Isend(&a[1*ROWS+2*pix],1,Column,W,0,new,&reqsend[2]);    // Send column west
	   MPI_Isend(&a[1*ROWS+JMAX],1,Column,E,0,new,&reqsend[3]); // Send column east
	   MPI_Isend(&a[1*ROWS+1*pix],1*pix,MPI_UNSIGNED_CHAR,NW,0,new,&reqsend[4]);
	   MPI_Isend(&a[1*ROWS+JMAX],1*pix,MPI_UNSIGNED_CHAR,NE,0,new,&reqsend[5]);
	   MPI_Isend(&a[IMAX*ROWS+1*pix],1*pix,MPI_UNSIGNED_CHAR,SW,0,new,&reqsend[6]);
	   MPI_Isend(&a[IMAX*ROWS+JMAX],1*pix,MPI_UNSIGNED_CHAR,SE,0,new,&reqsend[7]);


	   MPI_Irecv(&a[0*ROWS+1*pix], 1, Row, N, 0, new, &reqrecv[0]);  // recv from north
	   MPI_Irecv(&a[(IMAX+1)*ROWS+1*pix], 1, Row, S, 0, new, &reqrecv[1]); // recv from south
	   MPI_Irecv(&a[1*ROWS+0], 1, Column, W, 0, new, &reqrecv[2]); // recv from west
	   MPI_Irecv(&a[1*ROWS+JMAX+1*pix], 1, Column, E, 0, new, &reqrecv[3]);  // recv from east
	   MPI_Irecv(&a[0*ROWS+0], 1*pix, MPI_UNSIGNED_CHAR, NW, 0, new, &reqrecv[4]);
	   MPI_Irecv(&a[0*ROWS+JMAX+1*pix], 1*pix, MPI_UNSIGNED_CHAR, NE, 0, new, &reqrecv[5]);
	   MPI_Irecv(&a[(IMAX+1)*ROWS+0], 1*pix, MPI_UNSIGNED_CHAR, SW, 0, new, &reqrecv[6]);
	   MPI_Irecv(&a[(IMAX+1)*ROWS+JMAX+1*pix], 1*pix, MPI_UNSIGNED_CHAR, SE, 0, new, &reqrecv[7]);

	   // calculate inner subarray
	   unsigned char** h;
	   h = malloc(3*sizeof(char*));
	   for(int i=0;i<3;i++){
	     h[i] = malloc(3*sizeof(unsigned char));
	     for(int j=0;j<3;j++){
	       h[i][j] = i;
	     }
	   }

	//   #pragma omp parallel
	 //  {
	 //  #pragma omp for collapse(2)
	   for(int i = 2; i < IMAX; i++)
	      for(int j = 2; j < JMAX; j++)
	          b[i*ROWS+j] = output(a,i,j,h,1,ROWS);

	 //   #pragma omp single
	    MPI_Waitall(8,reqrecv,MPI_STATUSES_IGNORE);

	  //  #pragma omp barrier

	 //   #pragma omp parallel for
	    for(int j = 1; j<JMAX; j++){
	      b[1*ROWS+j] = output(a,1,j,h,1,ROWS);  // first row
	      b[IMAX*ROWS+j] = output(a,IMAX,j,h,1,ROWS); // last row
	    }

	    /* now we start from 2 and end up to IMAX-1 because the first and
	    the last element of the first column is also calculated in the
	    first and last row */
	 //   #pragma omp parallel for
	    for(int i = 2; i<IMAX; i++){
	      b[i*ROWS+1*pix] = output(a,i,1,h,1,ROWS);        // first column
	      b[i*ROWS+JMAX] = output(a,i,JMAX,h,1,ROWS);  // last column
	    }
	//  }

	    MPI_Waitall(8,reqsend,MPI_STATUSES_IGNORE);
	   // if(equals(a,b,IMAX,JMAX)) break;
      break;
	    a = b;
	}

	// open and write output image
	MPI_File_open(new,filenameout,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&f);

	/* remove halo points from b */
	b = repeat(b,IMAX,JMAX,pix);
	offsetb = 0;
	offsetf = (my_rank/sqrt(comm_sz))*sqrt(comm_sz)*IMAX*JMAX + (my_rank%(int)sqrt(comm_sz))*JMAX;
	for(int i = 0; i < IMAX; i++){
	  MPI_File_write_at(f,offsetf,&b[offsetb],JMAX,MPI_UNSIGNED_CHAR,&status);
	  offsetf += JMAX*sqrt(comm_sz);
	  offsetb += JMAX;
	}
	MPI_File_close(&f);

    MPI_Barrier(new);
    double end = MPI_Wtime();
    if(my_rank==0) printf("Time elapsed for process %d is %f seconds.\n", my_rank, end-start);

   MPI_Finalize();
   return 0;
} /* main */
