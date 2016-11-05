#ifndef DUDESEQ2_HPP
#define DUDESEQ2_HPP

#include <string>
#include <vector>
#include <boost/unordered_map.hpp>
using namespace std;

const int int_size = 10;

double pi_mat[int_size][int_size];
double pi_arr[int_size*int_size];
double inv_pi_arr[int_size*int_size];

int loss_mat[int_size][int_size] = {
	{0, 1, 1, 1, 1, 1, 1, 1, 1, 1},
	{1, 0, 1, 1, 1, 1, 1, 1, 1, 1},
	{1, 1, 0, 1, 1, 1, 1, 1, 1, 1},
	{1, 1, 1, 0, 1, 1, 1, 1, 1, 1},
	{1, 1, 1, 1, 0, 1, 1, 1, 1, 1},
	{1, 1, 1, 1, 1, 0, 1, 1, 1, 1},
	{1, 1, 1, 1, 1, 1, 0, 1, 1, 1},
	{1, 1, 1, 1, 1, 1, 1, 0, 1, 1},
	{1, 1, 1, 1, 1, 1, 1, 1, 0, 1},
	{1, 1, 1, 1, 1, 1, 1, 1, 1, 0},
};

typedef struct
{
	int i_flag, o_flag, k_flag, pi_flag, n_flag;
	char *i_file, *o_file, *pi_file;
	int k, n;
} option_t;

typedef struct
{
	int **t;
	int **c;
	int **a;
	int **g;
} count_t;

typedef boost::unordered_map<string, int> ma_map_t;
typedef boost::unordered_multimap<string, int> sl_map_t;


void GetArgument(int argc, char **argv, option_t &user_arg);
void SetProbability(char *pi_file);
double StringToDouble(const char *str_ptr);
void Quantize(const string flow_info, vector<string> &z_int, vector< pair<int, int> > &int_pos);
void SetInput(int argc, char **argv, option_t &user_arg, 
			  vector<string> &z_id, vector<string> &z_int, vector< pair<int, int> > &int_pos, int &count, int &max_length, int &k);
void Usage();
void AddToSlave(const string z_int, const pair<int, int> &int_pos, 
				const int k, sl_map_t &m, string **&joint_list, const int seq_num);
void SetHeap(count_t &m_arr, const int unique_key, bool mode);
void AddToMaster(count_t &m_arr, sl_map_t &m);
string Denoise(const string z_int, const pair<int, int> &int_pos, 
			   const int k, count_t &m_arr, string **&joint_list, ma_map_t &context_map, const int seq_num);

extern "C" {
	// LU decomoposition of a general matrix
	void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);
	
	// generate inverse of a matrix given its LU decomposition
	void dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork, int *INFO);
}

void Inverse(double *arr, double *inv_arr, int N)
{
	int *IPIV = new int[N+1];
	int LWORK = N*N;
	double *WORK = new double[LWORK];
	int INFO;
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			inv_arr[i*N + j] = arr[i*N + j];
	}
	
	dgetrf_(&N, &N, inv_arr, &N, IPIV, &INFO);
	dgetri_(&N, inv_arr, &N, IPIV, WORK, &LWORK, &INFO);
	
	delete IPIV;
	delete WORK;
}

#endif
