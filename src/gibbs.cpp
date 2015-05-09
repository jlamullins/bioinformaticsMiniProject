#include "gibbs.h"

gibbs::gibbs()
{
	srand(time(NULL));
}

gibbs::~gibbs()
{
	delete [] positions;
}

// Main implementation
void gibbs::algorithm()
{

	int skip = rand() % num_sequences;

	construct_pwm(sequences, positions, motif_length, skip);

/*
	std::string consensus;

	std::string alphabet = "ACGT";

	for(int i = 0; i < PWM.getCols(); i++)
	{
		double max = 0;
		int max_letter = 0;
		for(int j = 0; j < 4; j++)
		{
			if(PWM[j][i] > max)
			{
				max = PWM[j][i];
				max_letter = j;
			}
		}

		consensus += alphabet[max_letter];
	}
*/

	double max_score = 0.0;
	int max_pos = 0;

	for(int i = 0; i < sequences[skip].size() - motif_length + 1; i++)
	{
		double currScore = 0;

		for(int j = 0; j < motif_length; j++)
		{
			if(sequences[skip][i + j] == 'A')
			{
				currScore += log(PWM[0][j] / 0.25);
			}
			else if(sequences[skip][i + j] == 'C')
			{
				currScore += log(PWM[1][j] / 0.25);
			}
			else if(sequences[skip][i + j] == 'G')
			{
				currScore += log(PWM[2][j] / 0.25);
			}
			else if(sequences[skip][i + j] == 'T')
			{
				currScore += log(PWM[3][j] / 0.25);
			}
		}

		if(currScore > max_score)
		{
			max_score = currScore;
			max_pos = i;
		}
	}

	positions[skip] = max_pos;
}

// Assumes equal size strings/sequences
// Currently this computes a PPM, a probabilyt matrix, rather than a weight matrix.
// It should function the same as a PWM, but we might need to alter it to get a PWM if
// Something goes wrong.
void gibbs::construct_pwm(std::vector<std::string>& sequences, int * positions, int motif_length, int skip_seq)
{
	matrixf PFM, PWM;

	int size = sequences[0].size();

	int num_sequences = sequences.size();

	// If we want to skip a sequence, num_sequences == total sequences - 1
	if(skip_seq >= 0 || skip_seq < sequences.size())
	{
		num_sequences -= 1;
	}

	PFM.resize(4, motif_length);

	PWM = PFM;

	int pos_ctr = 0; // Each sequence has a different starting position, located in positions array

	for(auto it = sequences.begin(); it != sequences.end(); it++)
	{
		if(pos_ctr == skip_seq)
		{
			pos_ctr++;
			continue;
		}

		for(int pos = 0; pos < motif_length; pos++)
		{
			if((*it)[positions[pos_ctr] + pos] == 'A')
				PFM[0][pos]++;
			else if((*it)[positions[pos_ctr] + pos] == 'C')
				PFM[1][pos]++;
			else if((*it)[positions[pos_ctr] + pos] == 'G')
				PFM[2][pos]++;
			else if((*it)[positions[pos_ctr] + pos] == 'T')
				PFM[3][pos]++;
		}

		pos_ctr++;
	}

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < motif_length; j++)
		{
			PWM[i][j] = PFM[i][j] / num_sequences;
		}
	}

	this->PWM = PWM;
}

void gibbs::get_sequences(std::string filename, std::string motiflength, std::vector<std::string> & sequences)
{
	std::fstream file;
	std::string temp;


	file.open(filename.c_str(), std::ios::in);

	file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	while(std::getline(file, temp))
	{
		sequences.push_back(temp);
	}

	num_sequences = sequences.size();

	positions = new int[num_sequences];

	file.close();

	std::fstream mlfile;
	mlfile.open(motiflength.c_str(), std::ios::in);

	char c;
	mlfile.get(c);

	motif_length = c - '0';

	mlfile.close();
}

void gibbs::write_motif_to_file(std::string filename, std::string motif_name, int length, matrixf matrix)
{
	std::fstream file;
	
	file.open(filename.c_str(), std::ios::out | std::ios::trunc);	
	
	file << ">" << motif_name << "\t" << length << "\n";

	for(int i = 0; i < length; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			file << matrix[i * 4 + j] << "\t";
		}

		file << "\n";
	}
	
	file << "<";

	file.close();	
}

matrixf::matrixf()
{
	rows = 0;
	cols = 0;
}

matrixf::matrixf(int row, int col)
{
	array = new double*[row];

	rows = row;
	cols = col;

	for(int i = 0; i < row; i++)
	{
		array[i] = new double[col];
		for(int j = 0; j < col; j++)
		{
			array[i][j] = 0;
		}
	}
}

matrixf::~matrixf()
{
	for(int i = 0; i < rows; i++)
	{
		delete [] array[i];
	}
	if(rows > 0)
		delete [] array;
}

// Resize doesn't store any data from the original matrix
void matrixf::resize(int row, int col)
{
	for(int i = 0; i < rows; i++)
	{
		delete [] array[i];
	}
	if(rows > 0)
		delete [] array;
	
	rows = row;
	cols = col;

	array = new double*[rows];

	for(int i = 0; i < rows; i++)
	{
		array[i] = new double[cols];
		for(int j = 0; j < cols; j++)
		{
			array[i][j] = 0;
		}
	}
}

void matrixf::print()
{
	for(int i = 0; i < this->rows; i++)
	{
		for(int j = 0; j < this->cols; j++)
		{
			std::cout << (*this).array[i][j] << " ";
		}

		std::cout << std::endl;
	}
}

int main(int argc, char** argv)
{
	if(argc < 3) 
	{
		std::cout << "Invalid. Run with ./gibbs sequences.fa motiflength.txt" << std::endl;
		return -1;
	}

	gibbs gib;

	std::string seq = argv[1];
	std::string ml = argv[2];

	gib.get_sequences(seq, ml, gib.sequences);

	int samples = 150;

	// Initial positions are chosen randomly
	int seq_len = gib.sequences[0].size();
	for(int i = 0; i < gib.num_sequences; i++)
	{
		gib.positions[i] = rand() % seq_len;
	}

	for(int i = 0; i < samples; i++)
	{
		gib.algorithm();

		for(int j = 0; j < gib.sequences.size(); j++)
		{
			std::cout << gib.positions[j] << " ";
		}

		std::cout << std::endl;
	}

	return 0;
}
