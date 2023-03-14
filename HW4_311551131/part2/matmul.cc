/**/
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
using namespace std;
// *********************************************
// ** ATTENTION: YOU CANNOT MODIFY THIS FILE. **
// *********************************************

// Read size of matrix_a and matrix_b (n, m, l) and whole data of matrixes from stdin
//
// n_ptr:     pointer to n
// m_ptr:     pointer to m
// l_ptr:     pointer to l
// a_mat_ptr: pointer to matrix a (a should be a continuous memory space for placing n * m elements of int)
// b_mat_ptr: pointer to matrix b (b should be a continuous memory space for placing m * l elements of int)
void construct_matrices(int *n_ptr, int *m_ptr, int *l_ptr,
                        int **a_mat_ptr, int **b_mat_ptr)
{
    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int input_array[3];
    // Get n, m, i for everyone.
    if(world_rank) {
        MPI_Bcast(input_array, 3, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        for (int i = 0; i < 3; i++) 
            cin >> input_array[i];
        MPI_Bcast(input_array, 3, MPI_INT, 0, MPI_COMM_WORLD);
    }
    *n_ptr = input_array[0];
    *m_ptr = input_array[1];
    *l_ptr = input_array[2];
    // malloc 2D arrayint acnt = m*i;
	int nm = (*n_ptr)*(*m_ptr);
    int ml = (*l_ptr)*(*m_ptr);
    *a_mat_ptr = (int*)calloc(nm, sizeof(int));
	*b_mat_ptr = (int*)calloc(ml, sizeof(int));
    if(world_rank == 0){
        for (int i = 0; i < nm; i++)
            cin >> *(*a_mat_ptr + i);
		for(int i = 0 ; i < *m_ptr ; i++){
			for(int j = 0 ; j < *l_ptr ; j++){
				cin >> (*(*b_mat_ptr+i+j*(*m_ptr)));
			}
		}
    }
	MPI_Bcast(*a_mat_ptr, nm, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(*b_mat_ptr, ml, MPI_INT, 0, MPI_COMM_WORLD);
}

void matrix_multiply(const int n, const int m, const int l,
                     const int *a_mat, const int *b_mat)
{
    
    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	// initialize begin end
    int begin = (n / world_size) * world_rank;
    int end = (n / world_size) * (world_rank + 1);
    if (world_rank == world_size-1) 
		end = n;
    // initialize C
    int c_size = n*l;
    int *cmatrix;
    cmatrix = (int*)calloc(c_size, sizeof(int));
    
    // compute
    int ccnt = begin*l;
    for (int i = begin; i < end; i++) {
        for (int j = 0; j < l; j++) {
			int acnt = m*i;
            int sum = 0;
            int bcnt = j*m;
            for (int k = 0; k < m; k++){
                sum += a_mat[acnt++] * b_mat[bcnt++];
            }
            cmatrix[ccnt++] = sum;
			}
    }

    // // final
    int *rtn;
    rtn = (int*)calloc(c_size, sizeof(int));
    MPI_Reduce(cmatrix, rtn, n*l, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    ccnt = 0;
    if(world_rank == 0){
		for(int i = 0 ; i < n ; i++){
			for(int j = 0 ; j < l ; j++){
				printf("%d ", rtn[ccnt++]);
			}
			putchar('\n');
		}
	}
}

// Remember to release your allocated memory
void destruct_matrices(int *a_mat, int *b_mat)
{
    int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	if(world_rank == 0){
		free(a_mat);
		free(b_mat);
	}
}