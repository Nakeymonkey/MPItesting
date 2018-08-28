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
/*comput the maximum number*/
float compute_max(float *array, int tasks);
float *initial(string textname);
/*find block significant values*/
float blockmaxscore(string blockname, int start, int end);
int chunksize;
int numtasks, taskid;
float *significantvalue = new float[numofcluster];
float *significantnodes = new float[numofcluster*numtasks];
float Rdata[9611008];
int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	float *Arr = initial("Glioma.arrays.txt");
	MPI_Scatter( Arr, chunksize, MPI_FLOAT, Rdata,
			  chunksize, MPI_FLOAT, 0, MPI_COMM_WORLD);	
		 significantvalue[0] = blockmaxscore("G2", 0, 7);
		 significantvalue[1] = blockmaxscore("G2 oligo", 7, 45);
		 significantvalue[2] = blockmaxscore("G3", 45, 64);
		 significantvalue[3] = blockmaxscore(" G3 oligo", 64, 76);
		 significantvalue[4] = blockmaxscore(" G4", 76, 153);
		 significantvalue[5] = blockmaxscore(" non-tumour", 153, 176);

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
float compute_max(float *array, int tasks) {
	float max1 = 0;
	for (int i = 0; i < tasks; i++) {
		if (array[i] > max1) {
			max1 = array[i];
		}
	}
	return max1;
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

float blockmaxscore(string blockname, int start, int end) {
	int f = 0;
	float** Bar = new float*[bars];
	for (int n = 0; n < bars; ++n) {
		Bar[n] = new float[chunksize / bars];
	}

	for (int t = 0; t < chunksize / bars; t++) {
		for (int n = 0; n < bars; n++) {
			Bar[n][t] = Rdata[f];
			f++;
		}
	}
	float blockscore = Bar[start][0];
	for (int t = 0; t < chunksize / bars; t++) {
		for (int n = start; n < end; n++) {
			if (Bar[n][t] > blockscore) {
				blockscore = Bar[n][t];
			}
		}
	}
	cout << blockname <<" : "<< blockscore<< "\t";
	return blockscore;
}
