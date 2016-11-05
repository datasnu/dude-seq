#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <sstream>
#include <gsl/gsl_blas.h>
#include <omp.h>
#include "DUDE-Seq-2.hpp"
#include <sys/time.h>

int main(int argc, char **argv)
{
	timeval start, finish, begin, end, elapsed;
	gettimeofday(&start, NULL);
	
	option_t user_arg;
	int k, num_seq, max_length, unique_key = 0;
	vector<string> z_id, z_int;
	vector< pair<int, int> > int_pos;
	string **joint_list, *x_hat_int;
	count_t m_arr;
	ma_map_t context_map;
	ofstream out_stream;
	
	SetInput(argc, argv, user_arg, z_id, z_int, int_pos, num_seq, max_length, k);
	
	joint_list = new string* [num_seq];
	for (int i = 0; i < num_seq; i++)
		joint_list[i] = new string [max_length];
	
	out_stream.open(user_arg.o_file);
	if ( out_stream.fail() )
	{
		cout << "Error: Cannot open output file." << endl;
		exit(EXIT_FAILURE);
	}
	
	
	gettimeofday(&begin, NULL);
	omp_set_num_threads(user_arg.n);
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
			AddToSlave(z_int[i], int_pos[i], k, sl_m, joint_list, i);
		}
		
		// Creating unique context list
		for (sl_map_t::iterator it = sl_m.begin(); it != sl_m.end(); it = ( sl_m.equal_range(it->first) ).second)
		{
			sl_context_map.insert( pair<string, int>(it->first, sl_unique_key) );
			sl_unique_key++;
		}
		
		// Creating unique count array
		SetHeap(sl_m_arr, sl_unique_key, true);
		AddToMaster(sl_m_arr, sl_m);
		
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
				
				for (int i = 0; i < int_size; i++)
				{
					m_arr.t[ma_idx][i] += sl_m_arr.t[sl_idx][i];
					m_arr.a[ma_idx][i] += sl_m_arr.a[sl_idx][i];
					m_arr.c[ma_idx][i] += sl_m_arr.c[sl_idx][i];
					m_arr.g[ma_idx][i] += sl_m_arr.g[sl_idx][i];
				}
			}
		}
		
		SetHeap(sl_m_arr, sl_unique_key, false);
	}
	gettimeofday(&end, NULL);	timersub(&end, &begin, &elapsed);
	cout << "Step 1\n" << "========== " << user_arg.n << " thread(s) num_seqing time: " << (elapsed.tv_sec) + (elapsed.tv_usec / 1e6) << "s" << "\n";
	
	
	// Denoising sequences
	gettimeofday(&begin, NULL);
	x_hat_int = new string [num_seq];
	
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < num_seq; i++)
	{
		x_hat_int[i] = Denoise(z_int[i], int_pos[i], k, m_arr, joint_list, context_map, i);
	}
	
	for (int i = 0; i < num_seq; i++)
	{
		string out_str;
		out_stream << z_id[i] << '\n' << int_pos[i].first << ' ' << int_pos[i].second;
		for (int j = 0; j < x_hat_int[i].length(); j++)
		{
			out_str += ' ';
			out_str += x_hat_int[i][j];
		}
		out_stream << out_str << '\n';
	}
	out_stream.close();
	gettimeofday(&end, NULL);	timersub(&end, &begin, &elapsed);
	cout << "Step 2\n" << "========== " << user_arg.n << " thread(s) denoising time: " << (elapsed.tv_sec) + (elapsed.tv_usec / 1e6) << "s" << "\n\n";
	
	
	// Deallocating memory
	SetHeap(m_arr, unique_key, false);
	delete [] x_hat_int;
	
	for (int i = 0; i < num_seq; i++)
		delete [] joint_list[i];
	delete [] joint_list;
	
	gettimeofday(&finish, NULL);	timersub(&finish, &start, &elapsed);
	cout << "Elapsed time: " << (elapsed.tv_sec) + (elapsed.tv_usec / 1e6) << "s" << "\n";
	
	exit(EXIT_SUCCESS);
}


string Denoise(const string z_int, const pair<int, int> &int_pos,
			   const int k, count_t &m_arr, string **&joint_list, ma_map_t &context_map, const int seq_num)
{
	string joint_mer;
	string x_hat_int = z_int;
	int context_idx, channel, noisy;
	char denoised;
	
	int begin = int_pos.first;
	int end = int_pos.second;
	int margin = ((int)(k / 4) + 1) * 4;
	
	/*for (int i = begin+margin; i < end-margin; i++)*/
	for (int i = begin+k; i < end-k; i++)
	{
		joint_mer = joint_list[seq_num][i];
		noisy = z_int[i] - '0';
		
		context_idx = ( context_map.find(joint_mer) )->second;
		
		// Matrix computation start
		double f_term[int_size] = {};
		double s_term[int_size*int_size] = {};
		double scores[int_size] = {};
		double curr_m[int_size];
		
		channel = i % 4;
		switch (channel)
		{
			case 0:
				for (int j = 0; j < int_size; j++)
					curr_m[j] = m_arr.t[context_idx][j];
				break;
			case 1:
				for (int j = 0; j < int_size; j++)
					curr_m[j] = m_arr.a[context_idx][j];
				break;
			case 2:
				for (int j = 0; j < int_size; j++)
					curr_m[j] = m_arr.c[context_idx][j];
				break;
			case 3:
				for (int j = 0; j < int_size; j++)
					curr_m[j] = m_arr.g[context_idx][j];
				break;
		}
		
		gsl_matrix_view M_V = gsl_matrix_view_array(curr_m, 1, int_size);
		gsl_matrix_view INV_PI_V = gsl_matrix_view_array(inv_pi_arr, int_size, int_size);
		gsl_matrix_view F_TERM_V = gsl_matrix_view_array(f_term, 1, int_size);
		
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &M_V.matrix, &INV_PI_V.matrix, 0.0, &F_TERM_V.matrix);
		
		for (int j = 0; j < int_size; j++)
		{
			for (int k = 0; k < int_size; k++)
				s_term[k*int_size + j] = loss_mat[k][j] * pi_mat[k][noisy];
		}
		gsl_matrix_view S_TERM_V = gsl_matrix_view_array(s_term, int_size, int_size);
		gsl_matrix_view SCORES_V = gsl_matrix_view_array(scores, 1, int_size);
		
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &F_TERM_V.matrix, &S_TERM_V.matrix, 0.0, &SCORES_V.matrix);
		// Matrix computation end
		
		int min_idx = 0, num_min = 0;
		double min_score = scores[0];
		for (int j = 0; j < int_size; j++)
		{
			if (scores[j] <= min_score)
			{
				min_idx = j;
				min_score = scores[j];
			}
		}
		for (int j = 0; j < int_size; j++)
		{
			if (scores[j] == min_score) num_min++;
			if (num_min > 1) { min_idx = noisy; break; }
		}
		
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
		cout << "Error: Not enough parameters." << endl;
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
	
	if (!user_arg.i_flag | !user_arg.o_flag | !user_arg.k_flag | !user_arg.pi_flag)
	{
		Usage();
		exit(EXIT_FAILURE);
	}
	else if (!user_arg.n_flag)
	{
// 		user_arg.n = omp_get_num_procs();
		user_arg.n = 1;
	}
	
	if ( user_arg.k < 1 | user_arg.k > 9 )
	{
		cout << "Error: 'k' have to be in range 0 < k < 9." << endl;
		exit(EXIT_FAILURE);
	}
}


void Usage()
{
	cout << "USAGE: ./DUDE-Seq-2 <ARGUMENTS>" << endl
	<< "-i <file name>\t: input file [FLOWGRAM]" << endl
	<< "-o <file name>\t: output file" << endl
	<< "-k <integer>\t: size of the window (2 recommended)" << endl
	<< "-p <file name>\t: [optional] noise mechanism file (*.pi)" << endl << endl
	<< "EXAMPLE" << endl << "./DUDE-Seq-2 -i input.int -o output.int -k 2 -p Titanium.pi" << endl;
	exit(EXIT_FAILURE);
}


void SetProbability(char *pi_file)
{
	ifstream in_stream;
	string pi_info, token;
	double value;
	int i = 0, j = 0;
	
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
	
	for (i = 0; i < int_size; i++)
	{
		for (j = 0; j < int_size; j++)
			pi_arr[i*int_size + j] = pi_mat[i][j];
	}
	Inverse(pi_arr, inv_pi_arr, int_size);
}


double StringToDouble(const char *str_ptr)
{
	bool negative = false;
	double value = 0.0;
	
	if (*str_ptr == '-')
	{
		negative = true;
		++str_ptr;
	}
	
	while (*str_ptr >= '0' && *str_ptr <= '9')
	{
		value = (value * 10.0) + (*str_ptr - '0');
		++str_ptr;
	}
	
	if (*str_ptr == '.')
	{
		double frac_part = 0.0;
		int pow_part = 0;
		++str_ptr;
		
		while (*str_ptr >= '0' && *str_ptr <= '9')
		{
			frac_part = (frac_part * 10.0) + (*str_ptr - '0');
			++str_ptr;
			++pow_part;
		}
		
		value += frac_part / pow(10.0, pow_part);
	}
	
	if (negative) value = -value;
	
	return value;
}


void Quantize(const string flow_info, vector<string> &z_int, vector< pair<int, int> > &int_pos)
{
	double begin, end, value;
	string this_int, token;
	
	istringstream iss(flow_info);
	
	// Storing position values
	getline(iss, token, ' ');
	istringstream(token) >> begin;
	getline(iss, token, ' ');
	istringstream(token) >> end;
	int_pos.push_back( pair<int, int>(begin, end) );
	
	// Quantizing intensity values
	while ( getline(iss, token, ' ') )
	{
		value = StringToDouble( token.c_str() );
		value = floor(value + 0.5);
		if (value > 9) value = 9;
		
		this_int += (char)('0' + (int)value);
	}
	
	z_int.push_back(this_int);
}


void SetInput(int argc, char **argv, option_t &user_arg, 
			  vector<string> &z_id, vector<string> &z_int, vector< pair<int, int> > &int_pos, int &num_seq, int &max_length, int &k)
{
	ifstream in_stream;
	string flow_info;
	num_seq = 0, max_length = 0;
	
	GetArgument(argc, argv, user_arg);
	SetProbability(user_arg.pi_file);
	
	in_stream.open(user_arg.i_file);
	if ( in_stream.fail() )
	{
		cout << "Error: Cannot open input file." << endl;
		exit(EXIT_FAILURE);
	}
	
	while ( getline(in_stream, flow_info) )
	{
		z_id.push_back(flow_info);
		
		getline(in_stream, flow_info);
		Quantize(flow_info, z_int, int_pos);
		
		if (z_int[num_seq].length() > max_length)
			max_length = z_int[num_seq].length();
		
		num_seq++;
	}
	in_stream.close();
	
	k = user_arg.k;
}


void AddToSlave(const string z_int, const pair<int, int> &int_pos, 
				const int k, sl_map_t &m, string **&joint_list, const int seq_num)
{
	string left_mer, right_mer, joint_mer, curr_int;
	int noisy, channel, offset;
	
	int begin = int_pos.first;
	int end = int_pos.second;
	int margin = ((int)(k / 4) + 1) * 4;
	
	/*for (int i = begin+margin; i < end-margin; i++)
	{
		channel = i % 4;
		if (channel == 0)
		{
			left_mer = z_int.substr(i-k, k);
			right_mer = z_int.substr(i+4, k);
			joint_mer = left_mer + right_mer;
		}
		noisy = z_int[i] - '0';
		
		offset = channel * 10 + 10;
		m.insert( pair<string, int>(joint_mer, noisy+offset) );
		joint_list[seq_num][i] = joint_mer;
	}*/
	
	for (int i = begin + k; i < end - k; i++)
	{
		channel = i % 4;
		
		left_mer = z_int.substr(i-k, k);
		right_mer = z_int.substr(i+1, k);
		joint_mer = left_mer + right_mer;
		noisy = z_int[i] - '0';
		
		offset = channel * 10 + 10;
		m.insert( pair<string, int>(joint_mer, noisy+offset) );
		joint_list[seq_num][i] = joint_mer;
	}
}


void SetHeap(count_t &m_arr, const int unique_key, bool mode)
{
	if (mode == true)
	{	
		m_arr.t = new int* [unique_key];
		m_arr.a = new int* [unique_key];
		m_arr.c = new int* [unique_key];
		m_arr.g = new int* [unique_key];
		
		for (int i = 0; i < unique_key; i++)
		{
			m_arr.t[i] = new int [int_size];
			m_arr.a[i] = new int [int_size];
			m_arr.c[i] = new int [int_size];
			m_arr.g[i] = new int [int_size];
			
			for (int j = 0; j < int_size; j++)
			{
				m_arr.t[i][j] = 0;
				m_arr.a[i][j] = 0;
				m_arr.c[i][j] = 0;
				m_arr.g[i][j] = 0;
			}
		}
	}
	else
	{
		for (int i = 0; i < unique_key; i++)
		{
			delete [] m_arr.t[i];
			delete [] m_arr.a[i];
			delete [] m_arr.c[i];
			delete [] m_arr.g[i];
		}
		
		delete [] m_arr.t;
		delete [] m_arr.a;
		delete [] m_arr.c;
		delete [] m_arr.g;
	}
}


void AddToMaster(count_t &m_arr, sl_map_t &m)
{
	int context = 0, noisy;
	pair<sl_map_t::iterator, sl_map_t::iterator> range;
	
	for (sl_map_t::iterator it = m.begin(); it != m.end(); it = ( m.equal_range(it->first) ).second)
	{
		range = m.equal_range( it->first );
		
		for (sl_map_t::iterator r_it = range.first; r_it != range.second; r_it++)
		{
			noisy = r_it->second;
			
			if (noisy >= 10 && noisy < 20)
				m_arr.t[context][noisy-10]++;
			else if (noisy >= 20 && noisy < 30)
				m_arr.a[context][noisy-20]++;
			else if (noisy >= 30 && noisy < 40)
				m_arr.c[context][noisy-30]++;
			else
				m_arr.g[context][noisy-40]++;
		}
		
		context++;
	}
}
