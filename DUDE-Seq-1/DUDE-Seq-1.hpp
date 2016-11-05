#ifndef DUDESEQ1_HPP
#define DUDESEQ1_HPP

#include <string>
#include <vector>
#include <boost/unordered_map.hpp>
using namespace std;

const int nt_size = 4;

const string nt_order = "ATGCN";
double pi_mat[nt_size][nt_size] = {
	{0.9895, 0.0016, 0.0017, 0.0072},
	{0.0021, 0.9919, 0.0030, 0.0029},
	{0.0017, 0.0018, 0.9944, 0.0021},
	{0.0042, 0.0021, 0.0017, 0.9920}
}; // Illumina
// double pi_mat[nt_size][nt_size] = {
// 	{0.8122, 0.0034, 0.0894, 0.0950},
// 	{0.0096, 0.8237, 0.0808, 0.0859},
// 	{0.1066, 0.0436, 0.7774, 0.0724},
// 	{0.0704, 0.0690, 0.0889, 0.7717}
// }; // Nanopore

double pi_arr[nt_size*nt_size];
double inv_pi_arr[nt_size*nt_size];

int loss_mat[nt_size][nt_size] = {
	{0, 1, 1, 1},
	{1, 0, 1, 1},
	{1, 1, 0, 1},
	{1, 1, 1, 0}
};

typedef struct
{
	int i_flag, o_flag, k_flag, pi_flag, n_flag;
	char *i_file, *o_file, *pi_file;
	int k, n;
} option_t;

typedef int** count_t;

typedef boost::unordered_map<string, int> ma_map_t;
typedef boost::unordered_multimap<string, int> sl_map_t;


void GetArgument(int argc, char **argv, option_t &user_arg);
void SetProbability(int pi_flag, char *pi_file);
string NtToInt(const string z_seq);
string IntToNt(const string z_int);
void SetInput(int argc, char **argv, option_t &user_arg, 
			  vector<string> &z_id, vector<string> &z_int, vector<string> &z_qual, int &num_seq, int &max_length, int &k);
void Usage();
void AddToLocalM(const string z_int, const int k, sl_map_t &m);
void SetHeap(count_t &m_arr, const int unique_key, bool mode);
void AddToGlobalM(count_t &m_arr, sl_map_t &m);
string Denoise(const string z_int, const string z_qual, const int k, count_t &m_arr, ma_map_t &context_map);


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
