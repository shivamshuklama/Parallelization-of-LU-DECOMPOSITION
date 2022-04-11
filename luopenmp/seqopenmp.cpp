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
			{
				lu_value += temp_row[k]*temp_column[k];
			}
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
	// getting the current time
	auto start_timer = high_resolution_clock::now();

	// checking for correct command line input
	if(argc < 3)
	{
		cout << "Please input size of matrix and number of threads\n";
		return 0;
	}

	sscanf(argv[1], "%d", &MATRIX_SIZE);
	sscanf(argv[2], "%d", &NUM_THREADS);

	// below format will help us parallelize initialization of matrices.
	// declaring matrices as two dimensional vector of pointer to doubles
	//vector<vector<double>* > mat(MATRIX_SIZE, new vector<double>(MATRIX_SIZE, 0));
	vector<vector<double>* > mat, upper, lower, mat_dup;
	// declaring permutation matrix as vector of pointers to integers
	vector<int> permutation(MATRIX_SIZE);

	struct drand48_data buffers; // Each thread has a buffer corresponding to it

	// intializing all four matrices
	// seed the random number generator
	srand48_r((long int)time(0), &buffers);	/// seed differently for each thread ###################
	for(int i=0; i<MATRIX_SIZE; i++)
	{
		mat.push_back(new vector<double>(MATRIX_SIZE, (0)));
		mat_dup.push_back(new vector<double>(MATRIX_SIZE, (0)));
	}
	for(int i=0; i<MATRIX_SIZE; i++)
		upper.push_back(new vector<double>(MATRIX_SIZE, (0)));
	for(int i=0; i<MATRIX_SIZE; i++)
	{
		lower.push_back(new vector<double>(MATRIX_SIZE, (0)));
		(*lower[i])[i] = (1);
	}
	for(int i=0; i<MATRIX_SIZE; i++)
		permutation[i] = i;
	
	for(int i=0; i<MATRIX_SIZE; i++)
	{
		for(int j=0; j<MATRIX_SIZE; j++)
		{
			drand48_r( &buffers, &((*mat[i])[j]) );
			(*mat_dup[i])[j] = (*mat[i])[j];
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
		
	
		vector<double>* temp_row_pointer = mat[k];
		mat[k] = mat[k_dash];
		mat[k_dash] = temp_row_pointer;
	
		auto start = (*lower[k]).begin();
		auto end = (*lower[k]).begin()+k;
		double temp_sub_vector_lower_matrix[k];
		copy(start, end, temp_sub_vector_lower_matrix);

		start = (*lower[k_dash]).begin();
		end = (*lower[k_dash]).begin() + k;
		double temp_sub_vector_lower_matrix_kdash[k-1];
		copy(start, end, temp_sub_vector_lower_matrix_kdash);

		double temp_sub_vector_upper_kth_row[MATRIX_SIZE - 1 - k];
		start = (*mat[k]).begin() + (k+1);
		end = (*mat[k]).end();
		double temp_sub_vector_mat_kth_row[MATRIX_SIZE - k - 1];
		copy(start, end, temp_sub_vector_mat_kth_row);
		
		int temp = permutation[k];
		permutation[k] = permutation[k_dash];
		permutation[k_dash] = temp;
		
		(*upper[k])[k] = max_element;
		
		for(int j=0; j<=k-1; j++)
		{
			(*lower[k])[j] = temp_sub_vector_lower_matrix_kdash[j];
			(*lower[k_dash])[j] = temp_sub_vector_lower_matrix[j];
		}
		
		for(int j=k+1; j<MATRIX_SIZE; j++)
		{
			(*lower[j])[k] = (*mat[j])[k]/max_element;
			(*upper[k])[j] = temp_sub_vector_mat_kth_row[j - (k+1)];
		}
		// coping kth row of upper matrix: constant in next for loop
		start = (*upper[k]).begin() + k+1;
		end = (*upper[k]).end();
		copy(start, end, temp_sub_vector_upper_kth_row);
		// finally updating the input matrix
		for(int i=k+1; i<MATRIX_SIZE; i++)
		{
			double operand1 = (*lower[i])[k];
			for(int j=k+1; j<MATRIX_SIZE; j++)
				(*mat[i])[j] -= ((operand1 * temp_sub_vector_upper_kth_row[j-(k+1)]));
		}
	}
	// getting the end time
	auto end_time = high_resolution_clock::now();
	// get the duration
	auto duration_time = duration_cast<microseconds>(end_time - start_timer);
	cout << "Time taken: " << duration_time.count()/1000000.0 << " seconds" << endl << endl;

	//print_error(mat_dup, permutation, lower, upper);
	end_time = high_resolution_clock::now();
	duration_time = duration_cast<microseconds>(end_time - start_timer);
	//cout << "Time taken: " << duration_time.count()/1000000.0 << " seconds" << endl;

	return 0;
}

