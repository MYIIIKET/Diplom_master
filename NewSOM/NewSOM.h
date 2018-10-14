#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip> 
#include <vector> 
#include <stack>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

typedef map<int, set<int>> CLUSTERS;
typedef map<int, int> CLUSTER;

class NewSOM
{
public:
	NewSOM(const int entersNumber, const int height, const int width, bool debug = false, const double min = -0.5, const double max = 0.5,
		const int topol = RECT, const int windowLength = 2000, const float eta0 = 0.1, const float etaF = 0.08, const float sigma0 = 0.6, const float sigmaF = 0.2,
		const float beta = 0.7, const int countWindows=7);
	NewSOM(string&);
	~NewSOM();
	void trainMap(const double*);
	int findBMU(const double*,const int&);
	void saveMap(const char*);
	void saveMap(const string&);
	void getNormalizeVector(const string&);
	void getNormalizeVector(const char*);
	void normalizeDatabyMinMax(double* data);
	void initWeights();
	const enum cell { RECT, HEXA };
	double getWidth() { return mapWidth; }
	double getHeight() { return mapHeight; }
	double** getUmatrix();
private:
	int windowsNumber;
	int countWindows;
	void printClusters();

	double** map;
	CLUSTERS clusters;
	int* neighborhood;
	int mapWidth;
	int mapHeight;
	int neuronNumber;
	int neighborsNumber;
	double diag;		//”¡–¿“‹, ≈—À» Õ≈ œŒ“–≈¡”≈“—ﬂ!!!!!!!!!!!!!!!!!!!!
	int variableNumber;
	short cellType;

	double* D;
	double* Distances;
	double**V;

	double* quantizationError;
	double* neuronUtility;
	int *lastUpdate;
	
	double bigSigma;
	int iteration;
	int T;
	float eta0;
	float etaF;
	float sigma0;
	float sigmaF;
	double* normalizeVectorMAX;
	double* normalizeVectorMIN;
	float beta;

	void deleteFirstUpdate();
	void computeSigmaEta(double&, double&, double&);
	void updateConnections(const int&, const int&,const double &, const double &);

	bool isDebug;

	bool isOrderState;
	int whichLastBigger;
	int countWhichBigger;

	void changeState();
	int* neuronCoord(const int&);
	double rectDist(const int&, const int&);
	double hexDist(const int&, const int&);
	double neuronDist(const int&, const int&);


	double neighborhoodFun(const int&,const int&,const double&);
	//double(*neighborhoodFun)(const double&, const double&, const double&);
	double computeNeuronUtility();
	double evklDist(const double*, const double *);

	double driftFunction();
	double lastDriftValue;

	void randomInitialization();
	void saveVD();
	int** getAdjacencyMatrix();
	int getClustersNumber();
	void getClusters();
	void findCluster(const int&, const int&,int**,int*);
	void fillStack(stack<int>&, const int&, int*,int*);
	void findSubClusters();
	bool subCluster(set<int>&, const int&, const int& );
	void findMaxDensity( vector<int>&, set<int>& );
	bool isDensityMax(const int&);
	void findNeighbors(const int&);
	void deleteClusters(vector<int>&);
	bool isClustersNeighbors(set<int>&, set<int>&);
	bool isNeuronsNeighbors(const int&, const int&);
	bool isMerge(double**, const int& I, const int&, set<int>&, set<int>&);
	void makeSubCluster(const int&, CLUSTER&);
	void titleNeigbors(const int&, CLUSTER&, stack<int>&);
};

