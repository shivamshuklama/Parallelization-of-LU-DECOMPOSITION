#include<iostream>
#include<omp.h>
#include<vector>
#include<chrono>
#include<time.h>
#include<algorithm>
#include<cmath>

using namespace std;
using namespace std::chrono;

int NUM_THREADS;
int MATRIX_SIZE;



void print_matrix(vector<vector<double>* > mat, int r, int c)
{
	for(int i=0; i<r; i++)
	{
		for(int j=0; j<c; j++)	cout << (*mat[i])[j] << " ";
		cout << endl;
	}
	cout << endl;
}
void print_vector(vector<int> permutation, int size)
{
	for(int i=0; i<size; i++)	cout << permutation[i] << " ";
	cout << endl << endl;
}


void print_error(vector<vector<double>*> &mat_dup, vector<int> &permutation, vector<vector<double>*> &lower, vector<vector<double>*> &upper)
{
	//Obtaining the permuted matrix
	vector<vector<double>*> permuted_mat;
	for(int i=0; i<MATRIX_SIZE; i++)
	{
		permuted_mat.push_back(mat_dup[permutation[i]]);
	}

	
	double error_value = 0;
	//Subtracting LU from A and finding the error value
	#pragma omp parallel for num_threads(4) schedule(static, 4)
	for(int j=0; j<MATRIX_SIZE; j++)
	{
		vector<double> temp_column, temp_row;
		double squared_sum =0;
		for(int i=0; i<MATRIX_SIZE; i++)
		{
			double temp_value = (*upper[i])[j];
			temp_column.push_back(temp_value);
		}
		for(int i=0; i<MATRIX_SIZE; i++)
		{
			temp_row = *lower[i];
			double lu_value = 0;
			for(int k=0; k<MATRIX_SIZE; k++)
				lu_value += temp_row[k]*temp_column[k];
			squared_sum += pow ( ((*permuted_mat[i])[j] - lu_value), 2.0);
		}
		#pragma omp critical
		error_value += sqrt(squared_sum);
	}

	printf ("Error is %lf\n",error_value);


}

void print_vector_double(vector<double> vec, int r1, int r2)
{
	for(int i=r1; i<=r2; i++)
		cout << vec[i] << " ";
	cout << endl << endl;
}

// This function finds the maximum element from a fix column and between raw index r1 nad r2
// If such non-zero element exists then it returns its (index, max_value)
// Else it returns (-1, 0).
pair<int, double> get_maximum_element_index(vector<vector<double>* > &mat, int c, int r1, int r2)
{
	double max_element=0, pos_max_element=0;
	int index=-1;
	for(int i= r1; i<= r2; i++)
	{
		if(abs((*mat[i])[c]) > max_element)
		{
			max_element = abs((*mat[i])[c]);
			pos_max_element = (*mat[i])[c];
			index = i;
		}
	}
	return make_pair(index, pos_max_element);
}







int main(int argc, char *argv[])
{
	
	auto start_timer = high_resolution_clock::now();

	if(argc < 3)
	{
		cout << "enter size of matrix and number of threads\n";
		return 0;
	}

	sscanf(argv[1], "%d", &MATRIX_SIZE);
	sscanf(argv[2], "%d", &NUM_THREADS);

	vector<vector<double>* > mat, upper, lower, mat_dup;
	vector<int> permutation(MATRIX_SIZE);

	struct drand48_data buffers[NUM_THREADS]; // Each thread has a buffer corresponding to it

	// intializing all four matrices
	#pragma omp parallel num_threads(NUM_THREADS)
	{
		// seed the random number generator
		srand48_r((long int)time(0)+(long int)omp_get_thread_num(), &buffers[omp_get_thread_num()%NUM_THREADS]);
		#pragma omp sections
		{
			#pragma omp section
			{
				for(int i=0; i<MATRIX_SIZE; i++)
				{
					mat.push_back(new vector<double>(MATRIX_SIZE, (0)));
					mat_dup.push_back(new vector<double>(MATRIX_SIZE, (0)));
				}
			}
			#pragma omp section
			{
				for(int i=0; i<MATRIX_SIZE; i++)
					upper.push_back(new vector<double>(MATRIX_SIZE, (0)));
			}
			#pragma omp section
			{
				for(int i=0; i<MATRIX_SIZE; i++)
				{
					lower.push_back(new vector<double>(MATRIX_SIZE, (0)));
					(*lower[i])[i] = (1);
				}
			}
			#pragma omp section
			{
				for(int i=0; i<MATRIX_SIZE; i++)
					permutation[i] = i;
			}
		}
		// Initializing mat: 1000000 elements with random doubleing point values
		#pragma omp for schedule(static, 4) collapse(2)
		for(int i=0; i<MATRIX_SIZE; i++)
		{
			for(int j=0; j<MATRIX_SIZE; j++)
			{
				#pragma omp critical										// data dependency
				{
					drand48_r( &buffers[omp_get_thread_num()%4], &((*mat[i])[j]) );
					(*mat_dup[i])[j] = (*mat[i])[j];
				}
			}
		}
	}

	// outermost loop: non parallelizable due to data dependecies across iterations
	for(int k=0; k< MATRIX_SIZE; k++)
	{
		pair<int, double> max_of_column;
		max_of_column = get_maximum_element_index(mat, k, k, MATRIX_SIZE-1);
		int k_dash = max_of_column.first;
		double max_element = max_of_column.second;
		// if singular matrix then stop
		if(k_dash == -1)	{ cout << "Singular matrix.\n"; return 0;}
		
		// swapping k and k_dash row of matrix mat
		// This operation can not be parallelized because second last for loop assumes
		// two rows are swapped. If we do it parallely then second last for loop may
		// read non-swapped elements.
		vector<double>* temp_row_pointer = mat[k];
		mat[k] = mat[k_dash];
		mat[k_dash] = temp_row_pointer;
		// coping first (k-1) elements from kth row of lower triangular matrix because
		// while swapping elements in the loop every time going to that position through
		// pointers and copy the content into some temporary variable is time consuming
		// instead copy the whole subvector content to a temp vector and use that in
		// swapping elements. This will reduce order(n^2) instruction.
		auto start = (*lower[k]).begin();
		auto end = (*lower[k]).begin()+k;
		double temp_sub_vector_lower_matrix[k];
		copy(start, end, temp_sub_vector_lower_matrix);

		start = (*lower[k_dash]).begin();
		end = (*lower[k_dash]).begin() + k;
		double temp_sub_vector_lower_matrix_kdash[k];
		copy(start, end, temp_sub_vector_lower_matrix_kdash);

		double temp_sub_vector_upper_kth_row[MATRIX_SIZE - 1 - k];
		start = (*mat[k]).begin() + (k+1);
		end = (*mat[k]).end();
		double temp_sub_vector_mat_kth_row[MATRIX_SIZE - k - 1];
		copy(start, end, temp_sub_vector_mat_kth_row);
		// starting a parallel section
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			// only master threads completes the following work
			#pragma omp master
			{
				// swapping k and k_dash elements of permutation vector
				int temp = permutation[k];
				permutation[k] = permutation[k_dash];
				permutation[k_dash] = temp;
				// updating upper triangular matrix's (k, k) entry
				(*upper[k])[k] = max_element;
			}
			
			
			// swapping lower triangular matrix elements
			#pragma omp for schedule(static, 4) nowait
			for(int j=0; j<=k-1; j++)
			{
				#pragma omp critical								// data dependency
				{
					(*lower[k])[j] = temp_sub_vector_lower_matrix_kdash[j];
					(*lower[k_dash])[j] = temp_sub_vector_lower_matrix[j];
				}
			}
			// updating fractions in lower triangular matrix
			// and elements of upper triangular matrix
			
			
			#pragma omp for schedule(static, 4)
			for(int j=k+1; j<MATRIX_SIZE; j++)
			{
				(*lower[j])[k] = (*mat[j])[k]/max_element;
				#pragma omp critical								// data dependency
				{
					(*upper[k])[j] = temp_sub_vector_mat_kth_row[j - (k+1)];
				}
			}
			
			// barrier
			#pragma omp single
			{
				// coping kth row of upper matrix: constant in next for loop
				auto start = (*upper[k]).begin() + k+1;
				auto end = (*upper[k]).end();
				copy(start, end, temp_sub_vector_upper_kth_row);
			}
			
			// barrier 
			// finally updating the input matrix
			#pragma omp for schedule(static, 4)
			for(int i=k+1; i<MATRIX_SIZE; i++)
			{
				double operand1 = (*lower[i])[k];
				for(int j=k+1; j<MATRIX_SIZE; j++)
					(*mat[i])[j] -= ((operand1 * temp_sub_vector_upper_kth_row[j-(k+1)]));
			}
		}
	}
	// getting the end time
	auto end_time = high_resolution_clock::now();
	// get the duration
	auto duration_time = duration_cast<microseconds>(end_time - start_timer);
	cout << "Time taken: " << duration_time.count()/1000000.0 << " seconds" << endl << endl << endl;
	//print_error(mat_dup, permutation, lower, upper);
	end_time = high_resolution_clock::now();
	duration_time = duration_cast<microseconds>(end_time - start_timer);
	//cout << "Time taken: " << duration_time.count()/1000000.0 << " seconds" << endl;
	
	
	//print_matrix(upper,MATRIX_SIZE,MATRIX_SIZE);
	//print_matrix(lower,MATRIX_SIZE,MATRIX_SIZE);
	//print_vector(permutation,MATRIX_SIZE);

	return 0;
}

