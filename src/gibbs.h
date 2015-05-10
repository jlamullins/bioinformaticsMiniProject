#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

#include <stdlib.h>
#include <time.h>
#include <math.h>

class matrixf
{
	public:
		matrixf();
		matrixf(int row, int col);
		~matrixf();

		void resize(int row, int col);

		matrixf operator*(const matrixf& other)
		{
			matrixf retmat(this->rows, other.cols);
			
			if(this->cols != other.rows)
			{
				std::cout << "Invalid dimensions.\n";
				return retmat;
			}

			double val = 0;


			for(int i = 0; i < this->rows; i++)
			{
				for(int j = 0; j < other.cols; j++)
				{
					val = 0;
					for(int k = 0; k < this->cols; k++)
					{
						val += array[i][k] * other.array[k][j];
					}

					retmat[i][j] = val;
				}
			}

			return retmat;
		};
		matrixf& operator=(const matrixf& other)
		{
			this->resize(other.rows, other.cols);
			for(int i = 0; i < other.rows; i++)
			{
				for(int j = 0; j < other.cols; j++)
				{
					array[i][j] = other.array[i][j];
				}
			}

			return *this;
		};
		matrixf operator+(const matrixf& other)
		{
			matrixf retmat(this->rows, this->cols);
			if(this->rows != other.rows || this->cols != other.cols)
			{
				std::cout << "Invalid dimensions.\n";
				return retmat;
			}

			for(int i = 0; i < this->rows; i++)
			{
				for(int j = 0; j < this->cols; j++)
				{
					retmat.array[i][j] = array[i][j] + other.array[i][j];
				}
			}

			return retmat;
		};
		double*	 operator[](int index){ return array[index]; };

		int getRows(){ return rows; };
		int getCols(){ return cols; };

		void print();

	private:
		double** array;
		int rows, cols;
};

class gibbs
{
	public:
		gibbs();
		~gibbs();

		void algorithm();

		double getCurrentScore();
		void get_sequences(std::string filename, std::string motiflength, std::vector<std::string> & sequences);
		void write_motif_to_file(std::string filename, std::string motif_name, int length, matrixf& matrix);
		void construct_pwm(std::vector<std::string>& sequences, int * positions, int motif_length, int skip_seq);
		
		std::vector<std::string> sequences;

		int motif_length;
		matrixf PWM, PFM;
		int * positions, num_sequences;
};
