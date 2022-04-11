#include<iostream>
#include<vector>
#include<chrono>

#include<algorithm>
#include<cmath>
#include<cstdlib>
#include<cmath>
#include<time.h>

#include <unistd.h>
#include <sys/syscall.h>
#define gettid() syscall(SYS_gettid)

using namespace std;
using namespace std::chrono;

const int NUM_THREADS = 1;
int MATRIX_SIZE;
const int CHUNK_SIZE = 4;
int shared_counter = 0;
void init(void *data);
void processing(void *data);



struct matdata
{
	double **mat, **lower, **upper, **mat_dup;
	int *permutation;
	matdata()
	{
		
		mat = 	  (double **)(new double[MATRIX_SIZE]);
		mat_dup = (double **)(new double[MATRIX_SIZE]);
		lower =   (double **)(new double[MATRIX_SIZE]);
		upper =   (double **)(new double[MATRIX_SIZE]);
		permutation = new int[MATRIX_SIZE];

		for(int i=0; i<MATRIX_SIZE; i++)
		{
			mat[i] = new double[MATRIX_SIZE];
			mat_dup[i] = new double[MATRIX_SIZE];
			lower[i] = new double[MATRIX_SIZE];
			upper[i] = new double[MATRIX_SIZE];
			lower[i][i] = 1;
		}
	}
};

struct modified
{
	matdata *data_pointer;
	int k;
	int k_dash;
	modified(matdata *data_pointer, int k, int k_dash)
	{
		this->data_pointer = data_pointer;
		this->k = k;
		this->k_dash = k_dash;

	}
};


int get_maximum_element_index(double** mat, int c, int r1, int r2)
{
	double max_element=0;
	int index=-1;
	for(int i= r1; i<= r2; i++)
	{
		double temp = abs(mat[i][c]);
		if(temp > max_element)
		{
			max_element = temp;
			index = i;
		}
	}
	if(max_element == 0)
	{
		cout << "singular matrix"<<endl;
		return -1;
	}
	return index;
}

void print_error(matdata* data_pointer)
{
	// retrieving important data structures
	double **mat_dup = data_pointer->mat_dup;
	double **upper = data_pointer->upper;
	double **lower = data_pointer->lower;
	int* permutation = data_pointer->permutation;
	for(int i=0; i<MATRIX_SIZE; i++)
	{
		delete data_pointer->mat[i];
	}
	double** permuted_mat = data_pointer->mat;
	for(int i=0; i<MATRIX_SIZE; i++)
	{
		permuted_mat[i] = (mat_dup[permutation[i]]);
	}

	double error_value = 0;
	//Subtracting LU from A and finding the error value
	for(int j=0; j<MATRIX_SIZE; j++)
	{
		double temp_column[MATRIX_SIZE], temp_row[MATRIX_SIZE];
		double squared_sum =0;
		for(int i=0; i<MATRIX_SIZE; i++)
		{
			double temp_value = upper[i][j];
			temp_column[i] =(temp_value);
		}
		for(int i=0; i<MATRIX_SIZE; i++)
		{
			//temp_row = *(lower[i]);
			double lu_value = 0;
			for(int k=0; k<MATRIX_SIZE; k++)
			{
				lu_value += lower[i][k]*temp_column[k];
			}
			squared_sum += pow ( (permuted_mat[i][j] - lu_value), 2.0);
		}
		error_value += sqrt(squared_sum);
	}
	printf("Error is %lf\n",error_value);
	
	return;
}

// checks if elements from (k-1, k) onwords from kth row are same from mat and upper or not?
bool upper_mat_compatibility(double **upper, double **mat, int k)
{
	for(int j=k+1; j<MATRIX_SIZE; j++)
		if(mat[k][j] != upper[k][j])	return false;
	return true;
}

void print_matrix(double **mat, int size)
{
	for(int i=0; i<size; i++)
	{
		for(int j=0; j<size; j++)	cout << mat[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}
void print_vector(int *vec, int size)
{
	for(int i=0; i<size; i++)	cout << vec[i] << " ";
	cout << endl;
	cout << endl;
}



int main(int argc, char *argv[])
{
	// getting the current time
	auto start_timer = high_resolution_clock::now();

	// checking for correct command line input
	if(argc < 2)
	{
		cout << "give size_of_matrix\n";
		return 0;
	}

	sscanf(argv[1], "%d", &MATRIX_SIZE);

	// matdata has all matrices within it.
	matdata data = matdata();
	matdata *data_pointer = &data;

	init((void *)data_pointer);

	for(int k=0; k<MATRIX_SIZE; k++)
	{
		int k_dash = get_maximum_element_index(data.mat, k, k, MATRIX_SIZE-1);
		if(k_dash == -1)
		{
			cout << "Singular matrix.\n";
			return 0;
		}
		//swap permutation values
		int temp_permutation = data.permutation[k];
		data.permutation[k] = data.permutation[k_dash];
		data.permutation[k_dash] = temp_permutation;
		//swap k and k_dash rows of matrix mat
		double* temp_row_pointer = data.mat[k];
		data.mat[k] = data.mat[k_dash];
		data.mat[k_dash] = temp_row_pointer;
		//assign upper[k][k] the maximum_element
		data.upper[k][k] = data.mat[k][k];

		//creating threads for parallel computations
		modified data_p = modified(&data, k, k_dash);
		modified *data_pointer_p = &data_p;
		processing((void *)data_pointer_p);
		//joining threads

		// testing
		if(not upper_mat_compatibility(data.upper, data.mat, k))
			cout << "iteration: " << k << endl;

		shared_counter = 0;
	}

	// getting the end time
	auto end_time = high_resolution_clock::now();
	// get the duration
	auto duration_time = duration_cast<microseconds>(end_time - start_timer);
	cout << "Time taken: " << duration_time.count()/1000000.0 << " seconds" << endl;

	//ERROR value
	//print_error(data_pointer);
	end_time = high_resolution_clock::now();
	duration_time = duration_cast<microseconds>(end_time - start_timer);
	//cout << "Time taken after error calculation: " << duration_time.count()/1000000.0 << " seconds" << endl;
	for(int i=0; i<MATRIX_SIZE; i++)
	{
		delete data.mat_dup[i];
		delete data.lower[i];
		delete data.upper[i];
	}
	delete data.mat;
	delete data.mat_dup;
	delete data.lower;
	delete data.upper;
	delete data.permutation;

	exit(0);
}

// implementing the following for intializing mat elements
// #pragma omp parallel for schedule(static) num_threads(NUM_THREADS)
// initialize mat matrix element
void init(void *arg)
{

	matdata *data_pointer = (matdata *)arg;
	double **mat = data_pointer->mat;
	double **mat_dup = data_pointer->mat_dup;
	int* permutation = data_pointer->permutation;

	int tid = 0;

	struct drand48_data buffer;
	srand48_r((long int)time(0)+(long int)tid, &buffer);
	int start_iter_index = tid * ceil( ( (MATRIX_SIZE*(1.0)) / NUM_THREADS) );
	int end_iter_index = start_iter_index + ceil( ( (MATRIX_SIZE*(1.0)) / NUM_THREADS) ) - 1;
	end_iter_index = end_iter_index >= MATRIX_SIZE ? (MATRIX_SIZE-1) : end_iter_index;

	for(int i=start_iter_index; i<= end_iter_index; i++)
	{
		permutation[i] = i;
		for(int j=0; j<MATRIX_SIZE; j++)
		{
			drand48_r( &buffer, &(mat[i][j]) );
			mat[i][j] *= 1000;
			mat_dup[i][j] = mat[i][j];
		}
	}
}

void processing(void *arg)
{
	// retrieving important data structures
	modified *data_pointer_p = (modified *)arg;
	double **mat = data_pointer_p->data_pointer->mat;
	double **lower = data_pointer_p->data_pointer->lower;
	double **upper = data_pointer_p->data_pointer->upper;
	int k = data_pointer_p->k;
	int k_dash = data_pointer_p->k_dash;

	// getting current thread's thread_id
	int tid = 0;
	//calculating start and end indices
	int ceil_1 = ceil( ( (k*(1.0)) / NUM_THREADS ) );
	int ceil_2 = ceil( ( ((MATRIX_SIZE-k-1)*(1.0)) / NUM_THREADS ) );
	int start_iter_index_1 = tid * ceil_1;
	int end_iter_index_1 = start_iter_index_1 + ceil_1;
	end_iter_index_1 = end_iter_index_1  > (k) ? (k) : end_iter_index_1;
	int start_iter_index_2 = (k+1) + tid * ceil_2;
	int end_iter_index_2 = start_iter_index_2 + ceil_2;
	end_iter_index_2 = end_iter_index_2  > (MATRIX_SIZE) ? (MATRIX_SIZE) : end_iter_index_2;
	//swap k-1 elements of row numbers k and k_dash in the lower matrix
	for(int i=start_iter_index_1; i<end_iter_index_1; i++)
	{
		double temp_lower = lower[k][i];
		lower[k][i] = lower[k_dash][i];
		lower[k_dash][i] = temp_lower;
	}
	double u_k_k = upper[k][k];
	for(int i=start_iter_index_2; i<end_iter_index_2; i++)
	{
		lower[i][k] = mat[i][k]/u_k_k;
		upper[k][i] = mat[k][i];
	}
	for(int i=start_iter_index_2; i<end_iter_index_2; i++)
	{
		double operand1 = lower[i][k];
		for(int j=k+1; j<MATRIX_SIZE; j++)
			mat[i][j] -= operand1*upper[k][j];
	}
}


