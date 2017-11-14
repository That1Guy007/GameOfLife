#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define BORN 1
#define DIES 0

/* The Life function */
double life(int matrix_size, int ntimes, MPI_Comm comm){

  int      rank, size ; //rank is the node which is working on the program, the size is the total nodes in the system
  int      next, prev ;
  int      i, j, k;
  int      mysize, sum ;
  int    **matrix, **temp, **addr ;
  double   slavetime, totaltime, starttime ;
//  int      my_offset;

  /* Determine size and my rank in communicator */
  MPI_Comm_size(comm, &size) ;
  MPI_Comm_rank(comm, &rank) ;

  /* Set neighbors */
  if (rank == 0)
    prev = MPI_PROC_NULL;
  else
    prev = rank-1;
  if (rank == size - 1)
    next = MPI_PROC_NULL;
  else
    next = rank+1;


  /* Determine my part of the matrix */
  mysize = matrix_size/size ; //+ ((rank < (matrix_size % size)) ? 1 : 0 ) ;
  /*
   * rows
   */
  /*my_offset = rank * (matrix_size/size);
  if (rank > (matrix_size % size))
	  my_offset += (matrix_size % size);
  else
	  my_offset += rank;*/

  /* allocate the memory dynamically for the matrix */
  matrix = (int **)malloc(sizeof(int *)*(mysize+2)) ;
  temp = (int **)malloc(sizeof(int *)*(mysize+2)) ;

  for (i = 0; i < mysize+2; i++) {
    matrix[i] = (int *)malloc(sizeof(int)*(matrix_size+2)) ;
    temp[i] = (int *)malloc(sizeof(int)*(matrix_size+2))  ; //columns
  }

  /* Initialize the boundaries of the life matrix *//*
  for (j = 0; j < matrix_size+2; j++)
    matrix[0][j] = matrix[mysize+1][j] = temp[0][j] = temp[mysize+1][j] = DIES ;
  for (i = 0; i < mysize+2; i++)
    matrix[i][0] = matrix[i][matrix_size+1] = temp[i][0] = temp[i][matrix_size+1] = DIES ;*/

  /* Initialize the life matrix */

  srand(time(NULL));

  for (i = 0; i < mysize; i++)  { //rows
    for (j = 0; j< matrix_size; j++){ //collumns
    	if ((rand() % 2) == 1)
    		matrix[i][j] = BORN ;
      else
        matrix[i][j] = DIES ;
    }
  }

  /* Play the game of life for given number of iterations */
  starttime = MPI_Wtime() ;

  for (k = 0; k < ntimes; k++) {

    MPI_Request      req[4];
    MPI_Status       status[4];

    /* Send and receive boundary information
    MPI_Isend(&matrix[1][0],matrix_size+2,MPI_INT,prev,0,comm,req);
    MPI_Irecv(&matrix[0][0],matrix_size+2,MPI_INT,prev,0,comm,req+1);
    MPI_Isend(&matrix[mysize][0],matrix_size+2,MPI_INT,next,0,comm,req+2);
    MPI_Irecv(&matrix[mysize+1][0],matrix_size+2,MPI_INT,next,0,comm,req+3);
    MPI_Waitall(4, req, status);*/
    for(int loop = 0; loop < matrix_size; loop ++ ){
    //MPI_Send(&matrix[loop][0], 1, MPI_INT, rank -1, 0, MPI_COMM_WORLD);//left side
    MPI_Send(&matrix[0][loop], 1, MPI_INT, rank -1, 0, MPI_COMM_WORLD);//top to array above
    MPI_Send(&matrix[0][loop], 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);//top to array below
    MPI_Send(&matrix[matrix_size -1][loop], 1, MPI_INT, rank -1, 0, MPI_COMM_WORLD);//bottom
    MPI_Send(&matrix[matrix_size -1][loop], 1, MPI_INT, rank +1, 0, MPI_COMM_WORLD);//bottom
    //MPI_Send(&matrix[loop][matrix_size -1], 1, MPI_INT, rank -1, 0, MPI_COMM_WORLD);//right

    //MPI_Recv();
    }
    /* For each element of the matrix ... */
    for (i = 0; i < mysize; i++) {
      for (j = 0; j < matrix_size; j++) {

        /* find out the value of the current cell */

    	  /*
    	   * Need to wrap!
    	   * only need to wrap on i = 0 and i = Matrix_size/size, and matrix_size/size -1, and
    	   * i = matrix size.
    	   *
    	   */
    	if(i == 0){
    		if(j == 0){
    			//corner info
    			//sum = matrix[][];
    		}
    		//send and recieve top border nieghbor info
    	}
    	else if(j == 0 || j == mysize -1){
    		//left and right side info
    	}
    	else if(j == matrix_size-1){
    		//bottom border info
    	}
    	else{
    		//inner matrix
    	}

    	  /*
        sum = matrix[i-1][j-1] + matrix[i-1][j] + matrix[i-1][j+1]
          + matrix[i][j-1] + matrix[i][j+1]
            + matrix[i+1][j-1] + matrix[i+1][j] + matrix[i+1][j+1] ;*/

        /* check if the cell dies or life is born */
        if (sum < 2 || sum > 3)
          temp[i][j] = DIES ;
        else if (sum == 3)
          temp[i][j] = BORN ;
        else
          temp[i][j] = matrix[i][j] ;
		printf("%d ", matrix[i][j]);
      }
      printf("\n");
    }
    printf("\n");

    /* Swap the matrices */
    addr = matrix ;
    matrix = temp ;
    temp = addr ;
  }

  /* Return the average time taken/processor */
  slavetime = MPI_Wtime() - starttime;
  MPI_Reduce (&slavetime, &totaltime, 1, MPI_DOUBLE, MPI_SUM, 0, comm);



  return (totaltime/(double)size);
}


int main(argc, argv)
int argc ;
char *argv[] ;
{
  int rank, N, iters ;
  double time ;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

  /* If I'm process 0, determine the matrix size and # of iterations */
  /* This relies on the MPI implementation properly flushing output
     that does not end in a newline.  MPI does not require this, though
     high-quality implementations will do this.
   */
  if ( rank == 0 ) {
    printf("Matrix Size : ") ;
    scanf("%d",&N) ;
    printf("Iterations : ") ;
    scanf("%d",&iters) ;
  }

  /* Broadcast the size and # of iterations to all processes */
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD) ;
  MPI_Bcast(&iters, 1, MPI_INT, 0, MPI_COMM_WORLD) ;

  /* Call the life routine */
  time = life ( N, iters, MPI_COMM_WORLD );

  /* Print the total time taken */
  if (rank == 0)
    printf("[%d] Life finished in %lf seconds\n",rank,time/100);

  //MPE_Close_graphics(&graph);
  MPI_Finalize();
}
