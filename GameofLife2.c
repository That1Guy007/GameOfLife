



#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define BORN 1
#define DIES 0

/* The Life function */
double life(int matrix_size, int ntimes, MPI_Comm comm) {

	int rank, size; //rank is the node which is working on the program, the size is the total nodes in the system
	int next, prev, procSect;
	int i, j = 0, k;
	int mysize, sum = 0;
	//int **matrix, **temp, **addr;
	int *matrix, *temp, *addr;
	/*
	 * Matrix is now LINEAR
	 * This will allow us to use MPI_Scatter
	 *
	 */


	double slavetime, totaltime, starttime;
    /*
     * int      my_offset;
     * don't need offset, i break off sections of array or grid by horizontal rectangles
     * so each node never has to care about an offset
     */

	/* Determine size and my rank in communicator */
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	/* Set neighbors */

	/*
	 * prev is the previous node or process
	 * next is the next node or process
	 */


	if (rank == 0)
		prev = size -1;
		//prev = MPI_PROC_NULL; // means there is no previous process
	else
		prev = rank - 1;
	if (rank == size - 1)
		next = 0;
		//next = MPI_PROC_NULL; // means there is no next process
	else
		next = rank + 1;

	/* Determine my part of the matrix */
	mysize = (matrix_size / size) ; // Interesting : + ((rank < (matrix_size % size)) ? 1 : 0 ) ;

	procSect = mysize * matrix_size; // amount being handled per node


	/*
	 * rows or sections of initial array that the process
	 * or the node cares about
	 */

	temp = malloc(sizeof(int) * procSect); //columns

	if(rank ==0) {
		matrix = malloc(sizeof(int) * matrix_size*matrix_size);

		addr = malloc(sizeof(int) * matrix_size*matrix_size);




	/* Initialize the life matrix */

	srand(time(NULL));

	for (i = 0; i < matrix_size * matrix_size; i++) { //rows
		//for (j = 0; j < matrix_size; j++) { //collumns
			if ((rand() % 2) == 1)
				matrix[i] = BORN;
			else
				matrix[i] = DIES;
			if(i % mysize == 0 && i != 0)
							printf("\n");
			printf("%d ", matrix[i]);
		//}

	}

	printf("\n");

	for (i = 0; i < 2; i++) {
		//for (j = 0; j < matrix_size; j++) {
			addr[i] = DIES;
			printf("%d ", addr[i]);
		//}
			//printf("\n");

	}

	printf("\n");

	}



	/* Play the game of life for given number of iterations */
	starttime = MPI_Wtime() ;
	int mod;
	int mSize = matrix_size - 1;
	int offsetP;
	int offsetN;
	int maxS = (matrix_size * matrix_size);
	int uneven  =0;

	for (k = 0; k < ntimes; k++) {


		MPI_Request      req[4]; // handles the requested data
		MPI_Status       status[4];

		if(matrix_size % size == 0){
		if(MPI_Scatter(matrix, procSect, MPI_INT, //send the amount of procSect of matrix
					temp, procSect, MPI_INT,		  // To temp buffer
					0, MPI_COMM_WORLD) != MPI_SUCCESS){
				perror("Scatter Error");
				exit(1);
			}
		}
		/*
		 * handles the odd number of nodes or odd amount of rows
		 */
		else{
			uneven = 1;
			if(MPI_Scatter(matrix, procSect, MPI_INT, //send the amount of procSect of matrix
								temp, procSect, MPI_INT,		  // To temp buffer
								0, MPI_COMM_WORLD) != MPI_SUCCESS){
							perror("Scatter Error");
							exit(1);
			}
			if(rank == 0){
			for(int s = maxS - matrix_size; s < maxS; s++ ){
				MPI_Isend(&matrix[s], 1, MPI_INT, size -1, 0 ,comm, req);
			}
			}
			if(rank == size -1){
				for(int s = procSect; s <procSect + matrix_size ; s++ ){
					MPI_Irecv(&temp[s], 1, MPI_INT, 0, 0 ,comm, req);
							}
			}
		}
		/*
		 * sending each node an extra row from the previous nodes section
		 * this will allow me to update each row that would normally be skipped.
		 * add onto addr which will be used to add up the cells
		 */
		for(int idc = 0; idc <= mSize; idc++){

			MPI_Isend(&temp[idc], 1, MPI_INT, next, 0 , comm, req+1);

			MPI_Irecv(&addr[idc], 1, MPI_INT, prev, 0, comm, req+ 2);

		}

		/* For each element of the matrix WHICH WILL BE WORKED ON ... the Inner portion, ignore border ... */
			for (j = matrix_size; j < procSect + matrix_size; j++) {
				 mod = j%matrix_size;
				 offsetP= ((j - mod)- matrix_size);
				 offsetN = ((j- mod) + matrix_size);
				 if(j < procSect){
				/*
				 *  mSize = matrix_size - 1; // the section of each row
				 *  mod = j%matrix_size;		// the index of the row
				 *  offsetP = (j - matrix_size);	//the previous row
				 *  offsetN = (j + matrix_size); //the next row... (minus j mod, + matrix size)
				 */

				/* find out the value of the current cell */
				 if(j >= (procSect - matrix_size)){
					 if(mod == 0 || mod == mSize){
						 /*
						  * bot, botR, right, topR, top, topL, left, botL
						  */
						 if(mod == 0){
				 						sum = (addr[mod] + addr[mod + 1]
											+ temp[j +1] + temp[offsetP + 1]
											+ temp[offsetP] + temp[offsetP + mSize]
											+ temp[j + mSize] + addr[mod + mSize]);
						 }
						 else{
							 /*
							  * bot, botR, right, topR, top, topL, left, botL
							  */
							 sum = (addr[mod] + addr[0]
								 + temp[j - mSize] + temp[offsetP]
								 + temp[offsetP + mSize] + temp[offsetP + mSize -1]
								 + temp[j-1] + addr[mSize -1]);
						 }
					 }
					 else{
						 /*
						  * bot, botR, right, topR, top, topL, left, botL
						  */
						 sum = (addr[mod] + addr[mod  + 1]
								+ temp[j +1] + temp[offsetP + mod + 1]
								+ temp[offsetP + mod] + temp[offsetP + mod -1]
								+ temp[j -1] + addr[mod -1]);
					 }
				 }

				 else if (mod == 0 || mod == mSize) {
						//corner info
						if (mod == 0) //beginning of the row
							sum = (temp[mSize+ offsetP] + temp[j + mSize]
									+ temp[mod +(offsetN)] + temp[mod + (offsetN + mSize)]
									+ temp[offsetP + mod] + temp[offsetP + mod + 1]
									+ temp[j + 1] + temp[offsetN + mod + 1]);
						/*
						 * top left, leftWrap + bottom + bot left + top + top right + right + bottom right
						 */
						else//end of row
							sum = (temp[offsetP + mod] + temp[j- mSize] + temp[offsetP + mod -mSize]
									+ temp[offsetN] + temp[offsetN + mod -1] + temp[offsetP + mod -1]
									+ temp[offsetN + mod] + temp[j -1]);
						// top + rightW + topR + botR
						// bottomL + top L + bot + left
					}

				else {
						//for the rest of row
						//top + top left + top right + left and right + bot + botL + botR
						sum = (temp[offsetP + mod] + temp[offsetP + mod -1] + temp[offsetP + mod + 1]
								+ temp[j -1] + temp[j +1]
								+ temp[offsetN + mod] + temp[offsetN + mod - 1]
								+ temp[offsetN + mod + 1]);

						/*
						 * This will not work for all cases, but for the assignment, 1000 by 1000 it will
						 * (1 row per node is out of bounds)
						 */
					}
			}
				 else{//the last rank if it has more info from uneven node or rows
					 if(rank == size -1 && uneven){

							/*
							 *  mSize = matrix_size - 1; // the section of each row
							 *  mod = j%matrix_size;		// the index of the row
							 *  offsetP = (j - matrix_size);	//the previous row
							 *  offsetN = (j + matrix_size); //the next row... (minus j mod, + matrix size)
							 */

							/* find out the value of the current cell */
								 if(mod == 0 || mod == mSize){
									 /*
									  * bot, botR, right, topR, top, topL, left, botL
									  */
									 if(mod == 0){
							 						sum = (addr[mod] + addr[mod + 1]
														+ temp[j +1] + temp[offsetP + 1]
														+ temp[offsetP] + temp[offsetP + mSize]
														+ temp[j + mSize] + addr[mod + mSize]);
									 }
									 else{
										 /*
										  * bot, botR, right, topR, top, topL, left, botL
										  */
										 sum = (addr[mod] + addr[0]
											 + temp[j - mSize] + temp[offsetP]
											 + temp[offsetP + mSize] + temp[offsetP + mSize -1]
											 + temp[j-1] + addr[mSize -1]);
									 }
								 }
								 else{
									 /*
									  * bot, botR, right, topR, top, topL, left, botL
									  */
									 sum = (addr[mod] + addr[mod  + 1]
											+ temp[j +1] + temp[offsetP + mod + 1]
											+ temp[offsetP + mod] + temp[offsetP + mod -1]
											+ temp[j -1] + addr[mod -1]);
								 }

					 }

					 }


				/* check if the cell dies or life is born */
					if (sum < 2 || sum > 3)
					matrix[j] = DIES;
				else if (sum == 3)
					matrix[j] = BORN;

				printf("%d|%d|%d ", matrix[j], sum, j - matrix_size);
				sum = 0;
			}
			printf("\n \n");

			}

	/* Return the average time taken/processor */
	slavetime = MPI_Wtime() - starttime;
	MPI_Reduce(&slavetime, &totaltime, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

	return (totaltime / (double) size);
}

int main(argc, argv)
	int argc;char *argv[]; {
	int rank, N, iters;
	double time;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* If I'm process 0, determine the matrix size and # of iterations */
	/* This relies on the MPI implementation properly flushing output
	 that does not end in a newline.  MPI does not require this, though
	 high-quality implementations will do this.
	 */
	if (rank == 0) {
		printf("Matrix Size : ");
		scanf("%d", &N);
		printf("Iterations : ");
		scanf("%d", &iters);
	}

	/* Broadcast the size and # of iterations to all processes */
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&iters, 1, MPI_INT, 0, MPI_COMM_WORLD);

	/* Call the life routine */
	time = life(N, iters, MPI_COMM_WORLD);

	/* Print the total time taken */
	if (rank == 0)
		printf("[%d] Life finished in %lf seconds\n", rank, time / 100);

	MPI_Finalize();
}
