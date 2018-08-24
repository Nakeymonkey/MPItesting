#include<iostream>
#include <string>
#include<fstream>
#include <assert.h>
#include <sstream>
#include <vector>
#include"mpi.h"
using namespace std;
vector<float> arr;

float compute_avg(float *array, int tasks) {
	float sum= 0.000000f;
	int i;
	for (i = 0; i < tasks; i++) {
		sum += array[i];
	}
	return sum / tasks;
}

int main(int argc, char *argv[])
{
	/*initial parameter*/
	MPI_Init(&argc, &argv);
	int numtasks, taskid,  chunksize;	
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	std::ifstream infile("Glioma.arrays.txt");	
	int w=0;
	if (infile.is_open()) {
		string line;
		string names;
		getline(infile, line);
		getline(infile, line);
		getline(infile, line);
		getline(infile, line);
		while (getline(infile, line)) {			
			istringstream ss(line.c_str());
			ss >> names;
			ss >> names;
			ss >> names;		
			float word=0.000000f ;
			while (ss >> word) {				
			arr.push_back(word);
			w++;
			}
		}
	}
	
	
	float *arr1 = (float *)malloc(sizeof(float) * w);
	assert(arr1 != NULL);
	for (int i = 0; i < w; i++) {
		arr1[i] = arr[i];
	}

		infile.close();

	     chunksize = w / numtasks;
		 float *Rdata = (float *)malloc(sizeof(float) * chunksize);
		 assert(Rdata != NULL);
		  MPI_Scatter( arr1, chunksize, MPI_FLOAT, Rdata,
			  chunksize, MPI_FLOAT, 0, MPI_COMM_WORLD);
		 
		  int f = 0;
		  int bars = 176;
		  float** Bar = new float*[bars];
		  for (int n = 0; n < bars;++n) {
			  Bar[n] = new float[chunksize/bars];					
			  }
		 
		  for (int t = 0; t < chunksize / bars; t++) {
			  for (int n = 0; n < bars; n++) {
				  Bar[n][t] = Rdata[f];
				  f++;
			  }
		  }

		  float G2score=Bar[0][0];
		  for (int t = 0; t < chunksize/bars; t++) {
			  for (int n = 0; n < 7; n++) {
				  if (Bar[n][t] > G2score) {
					  G2score = Bar[n][t];
				  }
			  }
		  }
		  
		  cout <<  G2score << "\t";
		  float G2oligoscore = Bar[7][0];
		  for (int t = 0; t < chunksize / bars; t++) {
			  for (int n = 7; n < 45; n++) {
				  if (Bar[n][t] > G2oligoscore) {
					  G2oligoscore = Bar[n][t];
				  }
			  }
		  }
		  cout << G2oligoscore << "\t";

		  float G3score = Bar[45][0];
		  for (int t = 0; t < chunksize / bars; t++) {
			  for (int n = 45; n < 64; n++) {
				  if (Bar[n][t] > G3score) {
					  G3score = Bar[n][t];
				  }
			  }
		  }
		  cout << G3score << "\t";

		  float G3oligoscore = Bar[64][0];
		  for (int t = 0; t < chunksize / bars; t++) {
			  for (int n = 64; n < 76; n++) {
				  if (Bar[n][t] > G3oligoscore) {
					  G3oligoscore = Bar[n][t];
				  }
			  }
		  }
		  cout << G3oligoscore << "\t";

		  float G4score = Bar[76][0];
		  for (int t = 0; t < chunksize / bars; t++) {
			  for (int n = 76; n < 153; n++) {
				  if (Bar[n][t] > G4score) {
					  G4score = Bar[n][t];
				  }
			  }
		  }
		  cout << G4score << "\t";

		  float nontumourscore = Bar[153][0];
		  for (int t = 0; t < chunksize / bars; t++) {
			  for (int n = 153; n < 176	; n++) {
				  if (Bar[n][t] > nontumourscore) {
					  nontumourscore = Bar[n][t];
				  }
			  }
		  }
		  cout << nontumourscore << "\n";

		  float *G2score_avgs = (float *)malloc(sizeof(float) * numtasks);
		  assert(G2score_avgs != NULL);

		  MPI_Gather(&G2score, 1, MPI_FLOAT, G2score_avgs, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		  float score1 = compute_avg(G2score_avgs, numtasks);

		  cout <<" G2score "<< score1 << "\n";



		  float *G2oligoscore_avgs = (float *)malloc(sizeof(float) * numtasks);
		  assert(G2oligoscore_avgs != NULL);

		  MPI_Gather(&G2oligoscore, 1, MPI_FLOAT, G2oligoscore_avgs, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		  float score2 = compute_avg(G2oligoscore_avgs, numtasks);

		  cout << " G2oligoscore " << score2 << "\n";



		  float *G3score_avgs = (float *)malloc(sizeof(float) * numtasks);
		  assert(G3score_avgs != NULL);

		  MPI_Gather(&G3score, 1, MPI_FLOAT, G3score_avgs, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		  float score3 = compute_avg(G3score_avgs, numtasks);

		  cout <<  " G3score " << score3 << "\n";



		  float *G3oligoscore_avgs = (float *)malloc(sizeof(float) * numtasks);
		  assert(G3oligoscore_avgs != NULL);

		  MPI_Gather(&G3oligoscore, 1, MPI_FLOAT, G3oligoscore_avgs, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		  float score4 = compute_avg(G3oligoscore_avgs, numtasks);

		  cout << " G3oligoscore " << score4 << "\n";


		  float *G4score_avgs = (float *)malloc(sizeof(float) * numtasks);
		  assert(G4score_avgs != NULL);

		  MPI_Gather(&G4score, 1, MPI_FLOAT, G4score_avgs, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		  float score5 = compute_avg(G4score_avgs, numtasks);

		  cout << " G4score " << score5 << "\n";



		  float *nontumourscore_avgs = (float *)malloc(sizeof(float) * numtasks);
		  assert(nontumourscore_avgs != NULL);

		  MPI_Gather(&nontumourscore, 1, MPI_FLOAT, nontumourscore_avgs, 1, MPI_FLOAT,0, MPI_COMM_WORLD);
		  float score = compute_avg(nontumourscore_avgs, numtasks);
		
		  cout << " nontumourscore " << score << "\n";
		 

		  cin.get();
		  cin.get();

		  free(Rdata);
		  free(Bar);


	
	MPI_Finalize();

	

	return 0;
}
