#include<iostream>
#include <string>
#include<fstream>
#include <assert.h>
#include <sstream>
#include <vector>
#include"mpi.h"
using namespace std;
vector<float> arr;
/*comput the maximum number*/
float compute_max(float *array, int tasks) {
	float max1=0;
	for (int i = 0; i < tasks; i++) {
	if (array[i] > max1) {
		max1 = array[i];
	}
}
	return max1;
}

int main(int argc, char *argv[])
{
	/*initial processor*/
	MPI_Init(&argc, &argv);
	int numtasks, taskid,  chunksize;	
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	/*Read file*/
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
	infile.close();
	/*convert vector arr into float *arr */
	float *arr1 = (float *)malloc(sizeof(float) * w);
	assert(arr1 != NULL);
	for (int i = 0; i < w; i++) {
		arr1[i] = arr[i];
	}
	     chunksize = w / numtasks;

	/*Scatter the datas into processors*/
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
/*get G2, G2 oligo, G3, G3 oligo, G4,non-tumour's maximum value in different blocks */
		  float G2score=Bar[0][0];
		  for (int t = 0; t < chunksize / bars; t++) {
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

/*Pick several significant nodes*/

		  float *G2score_max = (float *)malloc(sizeof(float) * numtasks);
		  assert(G2score_max != NULL);

		  MPI_Gather(&G2score, 1, MPI_FLOAT, G2score_max, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		  float score1 = compute_max(G2score_max, numtasks);

		  cout <<" G2score "<< score1 << "\n";



		  float *G2oligoscore_max = (float *)malloc(sizeof(float) * numtasks);
		  assert(G2oligoscore_max != NULL);

		  MPI_Gather(&G2oligoscore, 1, MPI_FLOAT, G2oligoscore_max, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		  float score2 = compute_max(G2oligoscore_max, numtasks);

		  cout << " G2oligoscore " << score2 << "\n";



		  float *G3score_max = (float *)malloc(sizeof(float) * numtasks);
		  assert(G3score_max != NULL);

		  MPI_Gather(&G3score, 1, MPI_FLOAT, G3score_max, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		  float score3 = compute_max(G3score_max, numtasks);

		  cout <<  " G3score " << score3 << "\n";



		  float *G3oligoscore_max = (float *)malloc(sizeof(float) * numtasks);
		  assert(G3oligoscore_max != NULL);

		  MPI_Gather(&G3oligoscore, 1, MPI_FLOAT, G3oligoscore_max, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		  float score4 = compute_max(G3oligoscore_max, numtasks);

		  cout << " G3oligoscore " << score4 << "\n";


		  float *G4score_max = (float *)malloc(sizeof(float) * numtasks);
		  assert(G4score_max != NULL);

		  MPI_Gather(&G4score, 1, MPI_FLOAT, G4score_max, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		  float score5 = compute_max(G4score_max, numtasks);

		  cout << " G4score " << score5 << "\n";



		  float *nontumourscore_max = (float *)malloc(sizeof(float) * numtasks);
		  assert(nontumourscore_max != NULL);

		  MPI_Gather(&nontumourscore, 1, MPI_FLOAT, nontumourscore_max, 1, MPI_FLOAT,0, MPI_COMM_WORLD);
		  float score = compute_max(nontumourscore_max, numtasks);
		
		  cout << " nontumourscore " << score << "\n";
		 

		  cin.get();
		  cin.get();
		  /*clean arrays*/
		  free(Rdata);
		  free(Bar);


	
	MPI_Finalize();

	

	return 0;
}
