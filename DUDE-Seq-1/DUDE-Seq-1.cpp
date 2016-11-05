#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <gsl/gsl_blas.h>
#include <omp.h>
#include "DUDE-Seq-1.hpp"
#include <sys/time.h>

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char **argv)
{
	timeval start, finish, begin, end, elapsed;
	gettimeofday(&start, NULL);
	
	option_t user_arg;
	int k, num_seq, max_length, unique_key = 0;
	vector<string> z_id, z_int, z_qual;
	string z_qual_non = "";
	string **joint_list, *x_hat_int;
	count_t m_arr;
	ma_map_t context_map;
	ofstream out_stream;
	
	SetInput(argc, argv, user_arg, z_id, z_int, z_qual, num_seq, max_length, k);
	
	out_stream.open(user_arg.o_file);
	if ( out_stream.fail() )
	{
		cout << "Error: Cannot open output file." << endl;
		exit(EXIT_FAILURE);
	}
	
	
	gettimeofday(&begin, NULL);
	omp_set_num_threads(user_arg.n);
	cout << "Getting noise distribution..." << endl;
	#pragma omp parallel
	{
		sl_map_t sl_m;
		ma_map_t sl_context_map;
		int sl_unique_key = 0;
		count_t sl_m_arr;
		
		// Creating count multimap
		#pragma omp for schedule(static)
		for (int i = 0; i < num_seq; i++)
		{
			if (z_int[i].length() < (2*k + 1)) continue;
			AddToLocalM(z_int[i], k, sl_m);
		}
		
		// Creating unique context list
		for (sl_map_t::iterator it = sl_m.begin(); it != sl_m.end(); it = ( sl_m.equal_range(it->first) ).second)
		{
			sl_context_map.insert( pair<string, int>(it->first, sl_unique_key) );
			sl_unique_key++;
		}
		
		// Creating unique count array
		SetHeap(sl_m_arr, sl_unique_key, true);
		AddToGlobalM(sl_m_arr, sl_m);
		
		#pragma omp critical(GB_UNIQUE_CONTEXT)
		{
			for (ma_map_t::iterator it = sl_context_map.begin(); it != sl_context_map.end(); it++)
			{
				if ( context_map.find(it->first) == context_map.end() )
				{
					context_map.insert( pair<string, int>(it->first, unique_key) );
					unique_key++;
				}
			}
		}
		#pragma omp barrier
		
		#pragma omp master
		{
			SetHeap(m_arr, unique_key, true);
		}
		#pragma omp barrier
		
		#pragma omp critical(GB_M_ARR)
		{
			int ma_idx, sl_idx;
			
			for (ma_map_t::iterator it = sl_context_map.begin(); it != sl_context_map.end(); it++)
			{
				ma_idx = ( context_map.find(it->first) )->second;
				sl_idx = it->second;
				
				for (int j = 0; j < nt_size; j++)
					m_arr[ma_idx][j] += sl_m_arr[sl_idx][j];
			}
		}
		
		SetHeap(sl_m_arr, sl_unique_key, false);
	}
	gettimeofday(&end, NULL);	timersub(&end, &begin, &elapsed);
// 	cout << "Step 1\n" << "========== " << user_arg.n << " thread(s) counting time: " << (elapsed.tv_sec) + (elapsed.tv_usec / 1e6) << "s" << "\n";
	
	
	// Denoising sequences
	gettimeofday(&begin, NULL);
	x_hat_int = new string [num_seq];
	
	cout << "Correcting errors..." << endl << endl;
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < num_seq; i++)
	{
		if (z_int[i].length() < (2*k + 1))
		{
			x_hat_int[i] = z_int[i];
			continue;
		}
		
		if (z_qual.size() == 0)
			x_hat_int[i] = Denoise(z_int[i], z_qual_non, k, m_arr, context_map);
		else
			x_hat_int[i] = Denoise(z_int[i], z_qual[i], k, m_arr, context_map);
	}
		
	for (int i = 0; i < num_seq; i++)
	{
		out_stream << '>' << z_id[i] << endl;
		out_stream << IntToNt(x_hat_int[i]) << endl;
	}
	out_stream.close();
	gettimeofday(&end, NULL);	timersub(&end, &begin, &elapsed);
// 	cout << "Step 2\n" << "========== " << user_arg.n << " thread(s) denoising time: " << (elapsed.tv_sec) + (elapsed.tv_usec / 1e6) << "s" << "\n\n";
	
	
	// Deallocating memory
	SetHeap(m_arr, unique_key, false);
	delete [] x_hat_int;
	
	gettimeofday(&finish, NULL);	timersub(&finish, &start, &elapsed);
	cout << "Elapsed time: " << (elapsed.tv_sec) + (elapsed.tv_usec / 1e6) << "s" << endl;
	
	exit(EXIT_SUCCESS);
}


string Denoise(const string z_int, const string z_qual, const int k, count_t &m_arr, ma_map_t &context_map)
{
	string x_hat_int = z_int;
	string left_mer, right_mer, joint_mer;
	int noisy, context_idx;
	char denoised;
	
	for (int i = k; i < z_int.length()-k; i++)
	{
		left_mer = z_int.substr(i-k, k);
		right_mer = z_int.substr(i+1, k);
		joint_mer = left_mer + right_mer;
		noisy = z_int[i] - '0';
		
		context_idx = ( context_map.find(joint_mer) )->second;
		
		// Matrix computation start
		double f_term[nt_size] = {};
		double s_term[nt_size*nt_size] = {};
		double scores[nt_size] = {};
		
		double curr_m[nt_size];
		for (int j = 0; j < nt_size; j++)
			curr_m[j] = m_arr[context_idx][j];
		
		gsl_matrix_view M_V = gsl_matrix_view_array(curr_m, 1, nt_size);
		gsl_matrix_view INV_PI_V = gsl_matrix_view_array(inv_pi_arr, nt_size, nt_size);
		gsl_matrix_view F_TERM_V = gsl_matrix_view_array(f_term, 1, nt_size);
		
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &M_V.matrix, &INV_PI_V.matrix, 0.0, &F_TERM_V.matrix);
		
		for (int j = 0; j < nt_size; j++)
		{
			for (int sl_idx = 0; sl_idx < nt_size; sl_idx++)
				s_term[sl_idx*nt_size + j] = loss_mat[sl_idx][j] * pi_mat[sl_idx][noisy];
		}
		gsl_matrix_view S_TERM_V = gsl_matrix_view_array(s_term, nt_size, nt_size);
		gsl_matrix_view SCORES_V = gsl_matrix_view_array(scores, 1, nt_size);
		
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &F_TERM_V.matrix, &S_TERM_V.matrix, 0.0, &SCORES_V.matrix);
		// Matrix computation end
		
		int min_idx = 0, num_min = 0;
		double min_score = scores[0];
		for (int j = 0; j < nt_size; j++)
		{
			if (scores[j] <= min_score)
			{
				min_idx = j;
				min_score = scores[j];
			}
		}
		for (int j = 0; j < nt_size; j++)
		{
			if (scores[j] == min_score) num_min++;
			if (num_min > 1) { min_idx = noisy; break; }
		}
// 		to reduce FPs
// 		if ( (z_qual[i] - 33) > 30 ) { min_idx = noisy; }
		
		denoised = '0' + min_idx;
		x_hat_int[i] = denoised;
	}
	
	return x_hat_int;
}


void GetArgument(int argc, char **argv, option_t &user_arg)
{
	user_arg.i_flag = 0;
	user_arg.o_flag = 0;
	user_arg.k_flag = 0;
	user_arg.pi_flag = 0;
	
	if (argc < 2)
	{
		Usage();
		exit(EXIT_FAILURE);
	}
	
	for (int i = 1; i < argc; i++)
	{
		if ( !strcmp(argv[i], "-i") )
		{
			user_arg.i_flag = 1;
			user_arg.i_file = argv[++i];
		}
		else if ( !strcmp(argv[i], "-o") )
		{
			user_arg.o_flag = 1;
			user_arg.o_file = argv[++i];
		}
		else if ( !strcmp(argv[i], "-k") )
		{
			user_arg.k_flag = 1;
			istringstream(argv[++i]) >> user_arg.k;
		}
		else if ( !strcmp(argv[i], "-p") )
		{
			user_arg.pi_flag = 1;
			user_arg.pi_file = argv[++i];
		}
		else if ( !strcmp(argv[i], "-n") )
		{
			user_arg.n_flag = 1;
			istringstream(argv[++i]) >> user_arg.n;
		}
	}
	
	if (!user_arg.i_flag | !user_arg.o_flag | !user_arg.k_flag)
	{
		Usage();
		exit(EXIT_FAILURE);
	}
	
	if (!user_arg.n_flag)
	{
// 		user_arg.n = omp_get_num_procs();
		user_arg.n = 1;
	}
}


void Usage()
{
	cout << "USAGE: ./DUDE-Seq-1 <ARGUMENTS>" << endl
	<< "-i <file name>\t: input file [FASTA/FASTQ]" << endl
	<< "-o <file name>\t: output file" << endl
	<< "-k <integer>\t: size of the window (5 or 8 is recommended)" << endl
	<< "-p <file name>\t: [optional] noise mechanism file (*.pi)" << endl << endl
	<< "EXAMPLE" << endl << "./DUDE-Seq-1 -i input.fasta -o output.fasta -k 5 -p Illumina.pi" << endl;
	exit(EXIT_FAILURE);
}


void SetProbability(int pi_flag, char *pi_file)
{
	ifstream in_stream;
	string pi_info, token;
	double value;
	int i = 0, j = 0;
	
	if (pi_flag)
	{
		in_stream.open(pi_file);
		if ( in_stream.fail() )
		{
			cout << "Error: Cannot open mechanism file." << endl;
			exit(EXIT_FAILURE);
		}
		
		while ( getline(in_stream, pi_info) )
		{
			istringstream iss(pi_info);
			
			while ( getline(iss, token, ',') )
			{
				istringstream(token) >> value;
				pi_mat[i][j++] = value;
			}
			
			i++; j = 0;
		}
		in_stream.close();
	}
	                          
	for (i = 0; i < nt_size; i++)
	{
		for (j = 0; j < nt_size; j++)
			pi_arr[i*nt_size + j] = pi_mat[i][j];
	}
	Inverse(pi_arr, inv_pi_arr, nt_size);
}


void SetInput(int argc, char **argv, option_t &user_arg, 
			  vector<string> &z_id, vector<string> &z_int, vector<string> &z_qual, int &num_seq, int &max_length, int &k)
{
	ifstream in_stream;
	string nt_info;
	num_seq = 0, max_length = 0;
	
	gzFile file_ptr;
	kseq_t *z_str;
	
	GetArgument(argc, argv, user_arg);
	SetProbability(user_arg.pi_flag, user_arg.pi_file);
	
	// FASTA and FASTQ reader
	file_ptr = gzopen(user_arg.i_file, "r");
	z_str = kseq_init(file_ptr);
	
	while ( ( kseq_read(z_str) ) >= 0 )
	{
		z_id.push_back(z_str->name.s);
		z_int.push_back( NtToInt(z_str->seq.s) );
		if (z_str->qual.l) z_qual.push_back(z_str->qual.s);
		
		if (z_int[num_seq].length() > max_length)
			max_length = z_int[num_seq].length();
		
		num_seq++;
	}
	kseq_destroy(z_str);
	gzclose(file_ptr);
	
	k = user_arg.k;
}


string NtToInt(const string z_seq)
{
	string z_int;
	char curr_int;
	
	for (int i = 0; i < z_seq.length(); i++)
	{
		curr_int = '0' + nt_order.find(z_seq[i]);
		z_int += curr_int;
	}
	
	return z_int;
}


string IntToNt(const string z_int)
{
	string z_seq;
	int idx;
	
	for (int i = 0; i < z_int.length(); i++)
	{
		idx = z_int[i] - '0';
		z_seq += nt_order[idx];
	}
	
	return z_seq;
}


void AddToLocalM(const string z_int, const int k, sl_map_t &m)
{
	string left_mer, right_mer, joint_mer;
	int noisy;
	
	for (int i = k; i < z_int.length()-k; i++)
	{
		left_mer = z_int.substr(i-k, k);
		right_mer = z_int.substr(i+1, k);
		joint_mer = left_mer + right_mer;
		noisy = z_int[i] - '0';
		
		m.insert( pair<string, int>(joint_mer, noisy) );
	}
}


void SetHeap(count_t &m_arr, const int unique_key, bool mode)
{
	if (mode == true)
	{
		m_arr = new int* [unique_key];
		for (int i = 0; i < unique_key; i++)
		{
			m_arr[i] = new int [nt_size];
			for (int j = 0; j < nt_size; j++)
				m_arr[i][j] = 0;
		}
	}
	else
	{
		for (int i = 0; i < unique_key; i++)
			delete [] m_arr[i];
		delete [] m_arr;
	}
}


void AddToGlobalM(count_t &m_arr, sl_map_t &m)
{
	int context = 0, noisy;
	pair<sl_map_t::iterator, sl_map_t::iterator> key_range;
	
	for (sl_map_t::iterator it = m.begin(); it != m.end(); it = ( m.equal_range(it->first) ).second)
	{
		key_range = m.equal_range(it->first);
		
		for (sl_map_t::iterator r_it = key_range.first; r_it != key_range.second; r_it++)
		{
			noisy = r_it->second;
			m_arr[context][noisy]++;
		}
		
		context++;
	}
}
