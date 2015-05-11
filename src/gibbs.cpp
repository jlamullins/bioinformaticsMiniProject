#include "gibbs.h"

gibbs::gibbs()
{
	srand(time(NULL));
}

gibbs::~gibbs()
{
	delete [] positions;
}

double gibbs::getCurrentScore()
{
	construct_pwm(sequences, positions, motif_length, -1);

	double score = 0;
	int pos_ctr = 0;
	for(auto it = sequences.begin(); it != sequences.end(); it++)
	{
		for(int pos = 0; pos < motif_length; pos++)
		{
			if((*it)[positions[pos_ctr] + pos] == 'A')
				score += log(PFM[0][pos] / 0.25);
			else if((*it)[positions[pos_ctr] + pos] == 'C')
				score += log(PFM[1][pos] / 0.25);
			else if((*it)[positions[pos_ctr] + pos] == 'G')
				score += log(PFM[2][pos] / 0.25);
			else if((*it)[positions[pos_ctr] + pos] == 'T')
				score += log(PFM[3][pos] / 0.25);
		}

		pos_ctr++;
	}
}

// Main implementation
void gibbs::algorithm(int skip)
{

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

	this->PFM = PFM;
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

void gibbs::write_motif_to_file(std::string filename, std::string motif_name, int length, matrixf& matrix)
{
	std::fstream file;
	
	file.open(filename.c_str(), std::ios::out | std::ios::trunc);	
	
	file << ">" << motif_name << "\t" << length << "\n";

	for(int i = 0; i < matrix.getRows(); i++)
	{
		for(int j = 0; j < matrix.getCols(); j++)
		{
			file << matrix[i][j] << "\t";
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
	if(argc < 6) 
	{
		std::cout << "Invalid. Run with ./gibbs sequences.fa motiflength.txt predictedsites.txt predictedmotif.txt datasetNumber" << std::endl;
		return -1;
	}

	gibbs gib;

	std::string seq = argv[1];
	std::string ml = argv[2];
	std::string outsites = argv[3];
	std::string outmotif = argv[4];
	std::string dataset = argv[5];

	std::cout << "Getting sequences and motif length...";
	gib.get_sequences(seq, ml, gib.sequences);
	std::cout << " Done." << std::endl;

	// Usually reaches a constant set of positions by 100. To be verified
	int samples = 125;

	// Initial positions are chosen randomly
	int seq_len = gib.sequences[0].size();
	for(int i = 0; i < gib.num_sequences; i++)
	{
		gib.positions[i] = rand() % seq_len;
	}

	// Score each run, and pick positions from best score
	int num_runs = 55000 - (55000 % gib.sequences.size());	// Closest number to 25000 that's a multiple of # sequences
	double bestScore = 0, currentScore = 0;
	int * bestPositions = new int[gib.sequences.size()];	// Store best positions, given the total score of positions

	double percent = 0;

	for(int runs = 0; runs < num_runs; runs++)
	{
		for(int i = 0; i < gib.num_sequences; i++)
		{
			gib.positions[i] = rand() % seq_len;
		}

		for(int i = 0; i < samples; i++)
		{
			gib.algorithm(i % gib.sequences.size());

/*
			for(int j = 0; j < gib.sequences.size(); j++)
			{
				std::cout << gib.positions[j] << " ";
			}

			std::cout << std::endl;
*/
		}

		currentScore = gib.getCurrentScore();

		if(currentScore > bestScore)
		{
			bestScore = currentScore;

			for(int l = 0; l < gib.sequences.size(); l++)
			{
				bestPositions[l] = gib.positions[l];
			}
		}

		percent++;

		std::cout << "\rRunning Gibbs Sampling algorithm... " << percent / num_runs * 100 << "%" << "         \b\b\b\b\b\b\b\b\b";
	}

	std::cout << "\rRunning Gibbs Sampling algorithm...            \b\b\b\b\b\b\b\b\b\b\bDone." << std::endl;

	// Quite a lot of for loops :(
	for(int i = 0; i < gib.sequences.size(); i++)
	{
		gib.positions[i] = bestPositions[i];
	}

	std::string motifname = "MOTIF" + dataset;

	// Need to construct the PFM with no skips, which is why the last parameter
	// is -1
	gib.construct_pwm(gib.sequences, gib.positions, gib.motif_length, -1);

	gib.write_motif_to_file(outmotif, motifname, gib.motif_length, gib.PFM);

	/*
	for(int j = 0; j < gib.sequences.size(); j++)
	{
		std::cout << gib.positions[j] << " ";
	}
	std::cout << std::endl;
	*/

	std::fstream sitefile;

	sitefile.open(outsites.c_str(), std::ios::out | std::ios::trunc);

	for(int j = 0; j < gib.sequences.size(); j++)
	{
		sitefile << gib.positions[j] << " ";
	}

	sitefile << "\n";

	sitefile.close();

	return 0;
}
