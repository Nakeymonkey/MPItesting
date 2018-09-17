#include<iostream>
#include <string>
#include<fstream>
#include <assert.h>
#include <sstream>
#include <vector>
#include"mpi.h"
#include <math.h>

#define bars 176
#define numofcluster 6
using namespace std;
vector<float> arr;
/*initialize environment*/
float *initial(string textname);
/*comput the average number*/
void rowavgscore(int start, int end, int platform, int rows);
float correlationCoefficient(float X[], float Y[]);
/*initialize data*/
void dataInitial();
int chunksize;
int numtasks, taskid;
float *significantvalue = new float[numofcluster];
float *significantnodes = new float[numofcluster*numtasks];
float Rdata[9611008];
float** Bar = new float*[bars];
float** Rev = new float*[numofcluster];





int main(int argc, char *argv[])
{
	float sum = 0;
	MPI_Init(&argc, &argv);
	float *Arr = initial("Glioma.arrays.txt");


	MPI_Scatter(Arr, chunksize, MPI_FLOAT, Rdata,
		chunksize, MPI_FLOAT, 0, MPI_COMM_WORLD);

	dataInitial();
	
	ofstream myfile;
	myfile.open("correlationCoefficient.txt");
	
	for (int i = 0; i <chunksize / bars; i++) {
		float X[] = {Rev[0][i],Rev[1][i], Rev[2][i], Rev[3][i],  Rev[4][i], Rev[5][i] };
		for (int n = i; n <chunksize / bars; n++) {

			float Y[] = { Rev[0][n+1], Rev[1][n + 1], Rev[2][n + 1],  Rev[3][n + 1], Rev[4][n + 1], Rev[5][n+1] };
			myfile << correlationCoefficient(X, Y) << "\t";			
		}	
		cout << "\n";
		
	}
	myfile.close();
	
	MPI_Gather(significantvalue, numofcluster, MPI_FLOAT, significantnodes, numofcluster, MPI_FLOAT, 0, MPI_COMM_WORLD);

	cin.get();
	cin.get();
	/*clean arrays*/
	free(Rdata);

	MPI_Finalize();
	/*for (int i = 0; i < numofcluster*numtasks; i++) {
	cout << significantnodes[i] << "\t";
	}
	*/
	return 0;
}

float *initial(string textname) {

	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	std::ifstream infile(textname);
	int w = 0;
	if (infile.is_open()) {
		string line;
		while (getline(infile, line)) {
			istringstream ss(line.c_str());
			float word = 0.0000f;
			while (ss >> word) {
				arr.push_back(word);
				w++;
			}

		}

	}
	infile.close();

	float *arr1 = new float[w];
	for (int i = 0; i < w; i++) {
		arr1[i] = arr[i];
	}
	chunksize = w / numtasks;
	return arr1;
}

void rowavgscore(int start, int end, int platform, int rows) {
	float sum = 0;
	for (int t = start; t < end; t++) {
		sum = sum + Bar[t][rows];
	}
	Rev[platform][rows] = sum / (end - start);

}
void dataInitial() {
	int f = 0;

	for (int n = 0; n < bars; ++n) {
		Bar[n] = new float[chunksize / bars];
	}

	for (int t = 0; t < chunksize / bars; t++) {
		for (int n = 0; n < bars; n++) {
			Bar[n][t] = Rdata[f];
			f++;
		}
	}

	for (int n = 0; n < numofcluster; ++n) {
		Rev[n] = new float[chunksize / bars];
	}
	for (int n = 0; n < chunksize / bars; n++) {
		
			rowavgscore(0, 7, 0, n);
			rowavgscore(7, 45, 1, n);
			rowavgscore(45, 64, 2, n);
			rowavgscore(64, 76, 3, n);
			rowavgscore(76, 153, 4, n);
			rowavgscore(153, 176, 5, n);
		
	}
}

float correlationCoefficient(float X[], float Y[])
{

	float sum_X = 0, sum_Y = 0, sum_XY = 0;
	float squareSum_X = 0, squareSum_Y = 0;

	for (int i = 0; i < numofcluster; i++)
	{
		// sum of elements of array X.
		sum_X = sum_X + X[i];

		// sum of elements of array Y.
		sum_Y = sum_Y + Y[i];

		// sum of X[i] * Y[i].
		sum_XY = sum_XY + X[i] * Y[i];

		// sum of square of array elements.
		squareSum_X = squareSum_X + X[i] * X[i];
		squareSum_Y = squareSum_Y + Y[i] * Y[i];
	}
	

	float corr = (numofcluster * sum_XY - sum_X * sum_Y)/ sqrt((numofcluster * squareSum_X - sum_X * sum_X)*(numofcluster * squareSum_Y - sum_Y * sum_Y));

	return corr;
}
