#include<iostream>
#include <string>
#include<fstream>
#include <assert.h>
#include <sstream>
#include <vector>
#include"mpi.h"
#define bars 176
#define numofcluster 6
using namespace std;
vector<float> arr;
/*initialize environment*/
float *initial(string textname);
/*comput the average number*/
void rowavgscore(int start, int end, int platform,int rows);
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


	MPI_Scatter( Arr, chunksize, MPI_FLOAT, Rdata,
			  chunksize, MPI_FLOAT, 0, MPI_COMM_WORLD);	

	dataInitial();

	
	for (int i = 0; i < chunksize / bars; i++) {
		for (int n = 0; n < numofcluster; n++) {
			cout << Rev[n][i] << "\t";
		}
		cout << "\n";
	}

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
	//cout << arr[1];
	/*convert vector arr into float *arr */
	float *arr1 = new float[w];
	for (int i = 0; i < w; i++) {
		arr1[i] = arr[i];
	}
	chunksize = w / numtasks;
	return arr1;
}

void rowavgscore( int start, int end, int platform,int rows) {
	float sum = 0;	
		for (int t = start; t < end; t++) {
			sum = sum + Bar[t][rows];
		}
		Rev[platform][rows] = sum / (end-start);
	
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
		for (int i = 0; i < numofcluster; i++) {
			rowavgscore(0, 7, i, n);
			rowavgscore(7, 45, i, n);
			rowavgscore(45, 64, i, n);
			rowavgscore(64, 76, i, n);
			rowavgscore(76, 153, i, n);
			rowavgscore(153, 176, i, n);
		}
	}
}
