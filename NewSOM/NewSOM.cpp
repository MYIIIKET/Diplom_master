#include "NewSOM.h"
ofstream debugs;
ofstream qeHes;
ofstream clFile;

template<class T>
void printDebugs(string text, T* a, int size)
{
	debugs << text << endl << '\t';
	for (int i = 0; i < size; i++)
		debugs << a[i] << " ";
	debugs << endl;
	debugs<< endl;
}

template<class T>
void setArray(T*arr, int size, T val)
{
	for (int i = 0; i < size; i++)
		arr[i] = val;
}

NewSOM::NewSOM(const int entersNumber, const int height, const int width, bool d, const double min, const double max,
	const int topol /*=RECT*/, const int windowLength /*= 2000*/, const float eta0 /*= 0.1*/,
	const float etaF /*= 0.08*/, const float sigma0/* = 0.6*/, const float sigmaF /*= 0.2*/,
	const float beta /*= 0.7*/,const int countWindows) :
	variableNumber(entersNumber), mapHeight(height), mapWidth(width), cellType(topol),
	T(windowLength), eta0(eta0), etaF(etaF), sigma0(sigma0), sigmaF(sigmaF), beta(beta), isDebug(d),windowsNumber(countWindows)
{
	neuronNumber = mapWidth*mapHeight;
	this->countWindows = 1;
	map = new double*[neuronNumber];
	map[0] = new double[variableNumber*neuronNumber];
	for (size_t i = 1; i<neuronNumber; i++)
		map[i] = map[i - 1] + variableNumber;

	quantizationError = new double[windowLength];
	neuronUtility = new double[windowLength];

	lastUpdate = new int[neuronNumber];
	setArray(lastUpdate, neuronNumber, (-1)*T);

	if (cellType == RECT)
	{
		diag = sqrt(pow(mapWidth - 1, 2) + pow(mapHeight - 1, 2));
		neighborsNumber=4;
	}
	else
	{
		diag = neuronDist(0, neuronNumber - 1);
		neighborsNumber = 6;
	}

	D = new double[neuronNumber];
	//V = new double*[neuronNumber];
	//V[0] == new double[neuronNumber*neighborsNumber];
	//for (size_t i = 1; i<neuronNumber; i++)
	//	V[i] = V[i - 1] + neighborsNumber;
	//for (size_t i = 0; i < neuronNumber; i++)
	//	setArray(V[i], neighborsNumber, 0.0);

	V = new double*[2 * mapHeight - 1];
	V[0] = new double[(2 * mapHeight - 1)* (2 * mapWidth - 1)];
	for (size_t i = 1; i<(2 * mapHeight - 1); i++)
		V[i] = V[i - 1] + (2 * mapWidth - 1);
	for (size_t i = 0; i < (2 * mapHeight - 1); i++)
		setArray(V[i], (2 * mapWidth - 1), 0.0);

	setArray(D, neuronNumber, 0.0);

	normalizeVectorMAX = new double[entersNumber];
	setArray(normalizeVectorMAX, variableNumber, max);

	normalizeVectorMIN = new double[entersNumber];
	setArray(normalizeVectorMIN, variableNumber, min);

	Distances = new double[neuronNumber];
	setArray(Distances, neuronNumber, 0.0);

	bigSigma = 0;
	for (int i = 0; i < variableNumber; i++)
		bigSigma += pow((normalizeVectorMAX[i] - normalizeVectorMIN[i])*sqrt(variableNumber), 2);
	bigSigma = sqrt(bigSigma);

	neighborhood = new int[neighborsNumber];

	whichLastBigger = -2;
	countWhichBigger = 0;
	iteration = 0;
	isOrderState = true;
	randomInitialization();
	if (isDebug)
		debugs.open("DEBUG.txt", ios_base::out);
	//qeHe.open("qeHe.txt", ios_base::out);
}

NewSOM::NewSOM(string& path)
{
	string cell, neighbor;
	int mapSize1, mapSize2;
	ifstream fin(path);
	fin >> variableNumber >> cell >> mapSize1 >> mapSize2 >> neighbor;
	if (cell == "rect")cellType = RECT;
	else cellType = HEXA;

	normalizeVectorMAX = new double[variableNumber];
	normalizeVectorMIN = new double[variableNumber];

	//varName.resize(variableNumber);
	//fin >> varName[0];
	string test;
	for (int i = -1; i < variableNumber; i++)
		//*fin >> varName[i];
		fin >> test;

	mapWidth = mapSize1;
	mapHeight = mapSize2;

	neuronNumber = mapHeight*mapWidth;
	//neuronsClass.resize(neuronNumber);


	map = new double*[neuronNumber];
	map[0] = new double[variableNumber*neuronNumber];
	for (size_t i = 1; i<neuronNumber; i++)
		map[i] = map[i - 1] + variableNumber;
	//map = Matrix<double>(neuronNumber, entersNumber, 0);
	//this->entersNumber = entersNumber;


	for (int i = 0; i < mapHeight; i++)
	{
		for (int j = 0; j < mapWidth; j++)
			for (int k = 0; k < variableNumber; k++)
				fin >> map[mapHeight*j + i][k];
	}
	//for (int i = 0; i < neuronNumber; i++)
	//	neuronsClass[i] = '-';

	//for (int i = 0; i < neuronNumber; i++)
	//{
	//	for (int j = 0; j < variableNumber; j++)
	//		cout << map[i][j] << ' ';
	//	cout << endl;
	//}

	fin.close();
}



NewSOM::~NewSOM()
{
	delete[]map[0];
	delete[] map;
	delete[]V[0];
	delete[]V;
	delete[]D;
	delete[] Distances;
	delete[] normalizeVectorMAX;
	delete[] normalizeVectorMIN;
	delete[] quantizationError;
	delete[] lastUpdate;
	delete[] neuronUtility;
	delete[]neighborhood;
	debugs.close();
	qeHes.close();
}

double NewSOM::evklDist(const double* mapWeights, const double* data)
{
	double sum = 0;
	for (int i = 0; i < variableNumber; i++)
		sum += pow(mapWeights[i] - data[i], 2);
	return sqrt(sum);
}

double NewSOM::computeNeuronUtility()
{
	int num = 0;
	for (int i = 0; i < neuronNumber; i++)
	{
		if (lastUpdate[i] != (-1)*T)
			num++;
	}

	return (double)num / neuronNumber;
}


int NewSOM::findBMU(const double* data,const int& lastBMU)
{

	int index = 0;
	double minDist = 1000000000000000;
	for (int i = 0; i < neuronNumber; i++)
	{
		double dist = evklDist(map[i], data);
		if (dist < minDist && i!=lastBMU)
		{
			minDist = dist;
			index = i;
		}

		Distances[i] = dist;
	}

	if (isDebug) {
		printDebugs("quantizationError", quantizationError, T);
		printDebugs("BMU", &index, 1);
	}
	if (lastBMU == -1)
	{
		quantizationError[isOrderState ? iteration /*- 1*/ : T - 1] = minDist / bigSigma;//???????
		qeHes << quantizationError[isOrderState ? iteration/* - 1*/ : T - 1] << ' ';// << ' '<< neuronUtility[isOrderState ? iteration - 1 : T - 1] << endl;
	}
	return index;
}

int* NewSOM::neuronCoord(const int& neuron)
{
	int*  ans = new int[2];
	ans[0] = neuron / mapHeight;//x координата i-го нейрона;
	ans[1] = neuron %mapHeight;//y координата i-го нейрона;
	return ans;
}

double NewSOM::rectDist(const int& neuron, const int& winingNeuron)
{
	//double** coord=neuronCoord(neuron,winingNeuron);
	//double ans= sqrt(pow(coord[0][0] - coord[1][0], 2) + pow(coord[0][1] - coord[1][1], 2));//(x0-x1)^2+(y0-y1)^2
	int * nCoord = neuronCoord(neuron);
	int* wCoord = neuronCoord(winingNeuron);
	//double ans = sqrt(pow(nCoord[0] - wCoord[0], 2) + pow(nCoord[1] - wCoord[1], 2));
	double ans = 0;

	ans= abs(nCoord[0] - wCoord[0]) + abs(nCoord[1] - wCoord[1]);
	//ans = sqrt(pow(nCoord[0] - wCoord[0], 2) + pow(nCoord[1] - wCoord[1], 2));
	delete[] nCoord;
	delete[] wCoord;
	return ans;
}

double NewSOM::hexDist(const int& neuron, const int& winingNeuron)
{
	//double**coord = neuronCoord(neuron, winingNeuron);
	int* nCoord = neuronCoord(neuron);
	double xN = nCoord[0], yN = nCoord[1];
	int* wCoord = neuronCoord(winingNeuron);
	double xW = wCoord[0], yW = wCoord[1];

	if ((int)yN % 2 != 0)
		xN += 0.5;
	if ((int)yW % 2 != 0)
		xW += 0.5;
	yN *= sqrt(0.75);
	yW *= sqrt(0.75);
	//double ans = sqrt(pow(xN - xW, 2) + pow(yN - yW, 2));//(x0-x1)^2+(y0-y1)^2
	double ans = 0;

	//ans= abs(xN - xW) + abs(yN - yW);
	ans = sqrt(pow(xN - xW, 2) + pow(yN - yW, 2));//(x0-x1)^2+(y0-y1)^2

	delete[] nCoord;
	delete[] wCoord;

	return ans;
}

double NewSOM::neuronDist(const int& neuron, const int& winingNeuron)
{
	if (cellType == HEXA)return hexDist(neuron, winingNeuron);
	else return rectDist(neuron, winingNeuron);
}

double NewSOM::neighborhoodFun(const int& neuron, const int& winingNeuron,const double& sigma)
{
	
		//return (1/sigma)*exp((-1)*pow(neuronDist(neuron, winingNeuron) / (sigma*diag), 2));
	
		return exp((-1)*pow(neuronDist(neuron, winingNeuron) / (sigma*diag), 2));
}

template <class T>
void moveArray(T* a, int size)
{
	int i = size - 1;
	T val = a[i];
	for (i; i>0; i--)
	{
		a[i] += val;
		val = a[i] - val;
		a[i] = a[i] - val;
	}
	a[0] = val;
	//a[size-1] = (-1)*size - 1;

}

int findMin(int* arr, const int&len , const int& T)
{
	int a = 0;
	for (int i = 0; i < len; i++)
	{
		if ((arr[i] < a)&&T!= arr[i])
			a = arr[i];
	}
	return a;
}

void NewSOM::deleteFirstUpdate()
{
	for (int i = 0; i < neuronNumber; i++)
	{
		if (lastUpdate[i] < 0 && lastUpdate[i] * (-1) <= iteration && lastUpdate[i] * (-1) < T)
		{
			lastUpdate[i] = -T;
			continue;
		}
		if (lastUpdate[i] == 0 && iteration == 0)
		{
			lastUpdate[i] = -T;
			continue;
		}
	} 
	int a = findMin(lastUpdate, neuronNumber, -T);

}

int afk = 0;

void NewSOM::trainMap(const double* data)
{
	double etaI = 0, sigmaI = 0, neighborhood = 0, r =0;
	if (isDebug)printDebugs("ITERATION:", &iteration, 1);
	if (!isOrderState)
	{
		moveArray<double>(quantizationError, T);
		moveArray<double>(neuronUtility, T);
		if (isDebug) {
			printDebugs("neuronUtility", neuronUtility, T);
			printDebugs("quantizationError", quantizationError, T);
		}
		deleteFirstUpdate();
	}
	int BMU = findBMU(data,-1),
		BMU2 = findBMU(data, BMU);
	computeSigmaEta(sigmaI, etaI,r);
	r = 1 / (1 + exp(-1 * (double)(iteration+(countWindows-1)*T / (windowsNumber*T))));
	//r = 1 / (1 + exp(-1));
	qeHes << etaI << ' ' << sigmaI << ' ';
	for (int i = 0; i < neuronNumber; i++)
	{
		neighborhood = neighborhoodFun(i, BMU,sigmaI);
		if (neighborhood > 0.01)
		{
			for (int j = 0; j < variableNumber; j++)
				map[i][j] += etaI*neighborhood*(data[j] - map[i][j]);
			lastUpdate[i] = iteration;
				D[i] += r*exp((-0.5) * pow( /*evklDist(map[i],data)*/(Distances[i] /*/sigmaI*/), 2));
		}
		
	}

	neuronUtility[isOrderState ? iteration/*-1*/ : T - 1] = computeNeuronUtility();
	qeHes << neuronUtility[isOrderState ? iteration/*-1*/ : T - 1]<<' ';
	if (isDebug)
	{

		printDebugs("lastUpdate", lastUpdate, neuronNumber);
		printDebugs("neuronUtility", neuronUtility, T);
		for (int i = 0; i < mapHeight; i++)
		{
			for (int j = 0; j < mapWidth; j++)
			{
				debugs << lastUpdate[mapHeight*j + i] << ' ';
			}
			debugs << endl;
		}
		debugs << endl;

	}

	if (isNeuronsNeighbors(BMU, BMU2))qeHes << 0<< endl;
	else qeHes << 1 << endl;

	updateConnections(BMU, BMU2, r,etaI);
	
	iteration++;
	if (isOrderState)
	{
		if (iteration == T)
		{
			//saveMap(string("map")+std::to_string(afk++)+string(".cod"));
			changeState();
			countWindows++;
		}
	}
	else
	{

		if (/*countWhichBigger >= T-1||*/ countWhichBigger >= T)
		{
			//saveMap(string("map") + std::to_string(afk++) + string(".cod"));
			changeState();
		}
		if (iteration == T)
		{
			//saveMap(string("map") + std::to_string(afk++) + string(".cod"));
			for (int i = 0; i < neuronNumber; i++)
			{
				if (lastUpdate[i] > 0)
					lastUpdate[i] *= (-1);
			}
			//deleteFirstUpdate();
			if (whichLastBigger == -2)whichLastBigger = -1;//за один цикл лернинг не было случаем, когда дрифт больше оригина
			if (isDebug)printDebugs("end iterations T lastUpdate", lastUpdate, neuronNumber);
			iteration = 0;
			if(countWindows!=windowsNumber)
				countWindows++;
			else {
				printClusters();
				countWindows = 1;
			}
				

		}
	}
 }

void NewSOM::printClusters()
{
	//saveMap("map.cod");
	getClusters();
	int clusterNumber = getClustersNumber();
	saveVD();
	findSubClusters();
	clusterNumber = getClustersNumber();
	clFile << clusterNumber << endl;
	saveVD();


	for (size_t i = 0; i < (2 * mapHeight - 1); i++)
		setArray(V[i], (2 * mapWidth - 1), 0.0);
	setArray(D, neuronNumber, 0.0);

	for (int i = 0; i < clusters.size(); i++)
		clusters[i].clear();
	clusters.clear();
}

void NewSOM::findNeighbors(const int& ind)
{
	setArray(neighborhood, neighborsNumber, -5); 
	int x = ind / mapHeight, y = ind%mapHeight, nX, nY, neuronN;

	nY = y - 1;
	if (nY >= 0)
	{
		neighborhood[0] = x*mapHeight + nY;
	}
	nX = x + 1;
	if (nX<mapWidth)
	{
		neighborhood[1] = nX*mapHeight + y;
	}
	nY = y + 1;
	if (nY<mapHeight)
	{
		neighborhood[2] = x*mapHeight + nY;
	}
	nX = x - 1;
	if (nX >= 0)
	{
		neighborhood[3] = nX*mapHeight + y;
	}
	if (cellType == HEXA)
	{
		if (y % 2 == 0)
		{

			nY = y + 1;
			nX = x - 1;
			if (nX >= 0 && nY<mapHeight)
			{
				neighborhood[4] = nX*mapHeight + nY;
			}
			nY = y - 1;
			nX = x - 1;
			if (nX >= 0 && nY >= 0)
			{
				neighborhood[5] = nX*mapHeight + nY;
			}
		}
		else
		{

			nY = y - 1;
			nX = x + 1;
			if (nY >= 0 && nX<mapWidth)
			{
				neighborhood[4] = nX*mapHeight + nY;
			}

			nY = y + 1;
			nX = x + 1;
			if (nY<mapHeight&&nX<mapWidth)
			{
				neighborhood[5] = nX*mapHeight + nY;
			}
		}
	}
}

void NewSOM::updateConnections(const int& BMU1, const int& BMU2, const double& r,
	const double& etaI)
{
	int x = BMU1 / mapHeight, y = BMU1%mapHeight, nY = 0, nX = 0,
		neighborsNumber = 0;
	int neigh = isNeuronsNeighbors(BMU1, BMU2) ? neighborsNumber - 1 : neighborsNumber;
	double val =r/**neigh*/*etaI;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	nY = y - 1;
	neighborsNumber = x*mapHeight + nY;
	if (nY >= 0 && BMU2 != neighborsNumber)
	{
		V[y * 2 - 1][2 * x] -= V[y * 2 - 1][2 * x]*val;
	}
	else if (nY >= 0 && BMU2 == neighborsNumber)
	{
		V[y * 2 - 1][2 * x] -= (V[y * 2 - 1][2 * x]-1) * val;
	}

	nX = x - 1;
	neighborsNumber = nX*mapHeight + y;
	if (nX >= 0 && BMU2 != neighborsNumber)
	{
		V[2 * y][2 * x - 1] -= val*V[2 * y][2 * x - 1];
	}
	else if (nX >= 0 && BMU2 == neighborsNumber)
	{
		V[y * 2][2 * x - 1] -= val*(V[y * 2][2 * x - 1]-1);
	}

	nY = y + 1;
	neighborsNumber = x*mapHeight + nY;
	if (nY < mapHeight && BMU2 != neighborsNumber)
	{
		V[2 * y + 1][2 * x] -= val*V[2 * y + 1][2 * x];
	}
	else if (nY < mapHeight && BMU2 == neighborsNumber)
	{
		V[2 * y + 1][2 * x] -= val*(V[2 * y + 1][2 * x]-1);
	}

	nX = x + 1;
	neighborsNumber = nX*mapHeight + y;
	if (nX <mapWidth && BMU2 != neighborsNumber)
	{
		V[2 * y][2 * x + 1] -= val*V[2 * y][2 * x + 1];
	}
	else if (nX <mapWidth && BMU2 == neighborsNumber)
	{
		V[2 * y][2 * x + 1] -= val*(V[2 * y][2 * x + 1]-1);
	}

	if (cellType == HEXA)
	{
		if (y % 2 == 0)
		{
			nY = y - 1;
			nX = x - 1;
			neighborsNumber = nX*mapHeight + nY;
			if (nX >= 0 && nY >= 0 && BMU2 != neighborsNumber)
			{
				V[2 * y - 1][2 * x - 1] -= val*V[2 * y - 1][2 * x - 1];
			}
			else if(nX >= 0 && nY >= 0 && BMU2 == neighborsNumber)
			{
				V[2 * y - 1][2 * x - 1] -= val*(V[2 * y - 1][2 * x - 1] - 1);
			}
			nY = y + 1;
			nX = x - 1;
			neighborsNumber = nX*mapHeight + nY;
			if (nX >= 0 && nY < mapHeight && BMU2 != neighborsNumber)
			{
				V[2 * y + 1][2 * x - 1] -= val*V[2 * y + 1][2 * x - 1];
			}else if (nX >= 0 && nY < mapHeight && BMU2 == neighborsNumber)
			{
				V[2 * y + 1][2 * x - 1] -= val*(V[2 * y + 1][2 * x - 1] -1);
			}
		}
		else
		{
			nY = y - 1;
			nX = x + 1;
			neighborsNumber = nX*mapHeight + nY;
			if (nX < mapWidth && nY >= 0 && BMU2 != neighborsNumber)
			{
				V[2 * y - 1][2 * x + 1] -= val*V[2 * y - 1][2 * x + 1];
			}
			else if (nX < mapWidth && nY >= 0 && BMU2 == neighborsNumber)
			{
				V[2 * y - 1][2 * x + 1] -= val*(V[2 * y - 1][2 * x + 1] - 1);
			}
			nY = y + 1;
			nX = x + 1;
			neighborsNumber = nX*mapHeight + nY;
			if (nX < mapWidth && nY < mapHeight&& BMU2 != neighborsNumber)
			{
				V[2 * y + 1][2 * x + 1] -= val*V[2 * y + 1][2 * x + 1];
			
			}else if (nX < mapWidth && nY < mapHeight&& BMU2 == neighborsNumber)
			{
				V[2 * y + 1][2 * x + 1] -= val*(V[2 * y + 1][2 * x + 1] - 1);
			}
		}

	}
}

void NewSOM::computeSigmaEta(double& sigmaI, double& etaI,double &r)
{
	double arg = 0;
	if (isOrderState)
	{
		arg = (double)iteration / (T - 1);
		sigmaI = sigma0*pow(sigmaF / sigma0, arg);
		etaI = eta0*pow(etaF / eta0, arg);
		r = 1 / (1 + exp((-1)*arg));
	}
	else
	{

		arg = driftFunction();
		//r= 1 / (1 + exp(-1));
		if (isDebug)
			printDebugs("driftFunction", &arg, 1);

		if (arg < lastDriftValue)
		{
			sigmaI = (iteration == 0) && (whichLastBigger == -2) ? sigmaF : (arg*sigmaF) / lastDriftValue;
			etaI = (iteration == 0) && (whichLastBigger == -2) ? etaF : (arg*etaF) / lastDriftValue;
			r = (iteration == 0) && (whichLastBigger == -2) ? 1 / (1 + exp(-1)) : (arg* (1 / (1 + exp(-1)))) / lastDriftValue;
			if (countWhichBigger)
			{
				countWhichBigger = 0;
				whichLastBigger = -3;
			}
		}
		else
		{
			sigmaI = sigmaF;
			etaI = etaF;
			r = 1 / (1 + exp(-1));
			if (whichLastBigger == ((iteration - 1 == -1/*0*/) ? T - 1 : iteration - 1))
			{
				whichLastBigger = iteration;
				countWhichBigger++;
			}
			else
			{
				whichLastBigger = iteration;
				countWhichBigger = 1;
			}
			//countWhichBigger++;
		}
	}
}

void NewSOM::changeState()
{
	if (isOrderState)
	{
		lastDriftValue = driftFunction();
		for (int i = 0; i < neuronNumber; i++)
		{
			if (lastUpdate[i] != -T)lastUpdate[i] *= -1;
			//if (lastUpdate[i] ==0)lastUpdate[i] =lastVal;
		}
		if (isDebug) {
			printDebugs("--------------------------LEARNING STATE--------------------------", &lastDriftValue, 1);
			printDebugs("lastUpdate", lastUpdate, neuronNumber);
		}

		iteration = 0;
		countWhichBigger = 0;
		whichLastBigger = -2;
		isOrderState = false;
	}
	else
	{
		if (isDebug)
			printDebugs("--------------------------ORDER STATE--------------------------", &lastDriftValue, 1);
		setArray(lastUpdate, neuronNumber, -T);
		countWindows = 1;


		iteration = 0;
		countWhichBigger = 0;
		whichLastBigger = -2;
		isOrderState = true;
	}
}

double NewSOM::driftFunction()
{
	double qe = 0;
	for (int i = 0; i <T; i++)
		qe += quantizationError[i];
	qe /= T;
	double he = 0;
	for (int i = 0; i < T; i++)
	{
		he += neuronUtility[i];
	}
	he /= T;
	return beta*qe + (1 - beta)*(1 - he);
}

void NewSOM::initWeights()
{
	bigSigma = 0;
	for (int i = 0; i < variableNumber; i++)
		bigSigma += pow((normalizeVectorMAX[i] - normalizeVectorMIN[i]), 2);
	bigSigma = sqrt(bigSigma)*sqrt(variableNumber);
	randomInitialization();
}

void NewSOM::normalizeDatabyMinMax(double* data)
{
	for (int i = 0; i < variableNumber; i++)
	{
		if(normalizeVectorMAX[i]!=0)
			data[i] = (data[i] - normalizeVectorMIN[i]) / (normalizeVectorMAX[i] - normalizeVectorMIN[i]);
		else 
			data[i] = 0;
	}
}

void NewSOM::randomInitialization()
{
	/*srand(time(NULL));*/
	for (int i = 0; i < neuronNumber; i++)
		for (int j = 0; j < variableNumber; j++)
			map[i][j] = (double)(rand()) / RAND_MAX*(normalizeVectorMAX[j] - normalizeVectorMIN[j]) + normalizeVectorMIN[j];
}

void NewSOM::saveMap(const string &path)
{
	saveMap(path.c_str());
}

void NewSOM::saveMap(const char* path)
{
	ofstream fout(path, ios_base::out);
	fout << variableNumber << " ";
	if (cellType == RECT)fout << "rect";
	else fout << "hexa";
	fout << " " << mapWidth << " " << mapHeight << " ";
	fout << "gaussian";
	fout << endl;

	fout << "#n ";

	for (int i = 0; i < variableNumber; i++)
		fout << '-' << " ";
	fout << endl;

	for (int i = 0; i < mapHeight; i++)
	{
		for (int j = 0; j < mapWidth; j++)
		{
			for (int k = 0; k < variableNumber; k++)
				fout << map[mapHeight*j + i][k] << ' ';
			fout << endl;
		}
	}

	fout.close();
}

void NewSOM::getNormalizeVector(const string& path)
{
	getNormalizeVector(path.c_str());
}


double** NewSOM::getUmatrix()
{
	double*** M = new double**[variableNumber];
	M[0] = new double*[mapHeight*variableNumber];
	M[0][0] = new double[neuronNumber*variableNumber];

	for (int i = 1; i < variableNumber; i++)
		M[i] = M[i - 1] + mapHeight;

	for (int j = 1; j < mapHeight*variableNumber; j++)
		M[0][j] = M[0][j - 1] + mapWidth;

	for (int i = 0; i < variableNumber; i++)
		for (int j = 0; j < mapHeight; j++)
			for (int k = 0; k<mapWidth; k++)
				M[i][j][k] = map[k*mapHeight + j][i];

	int HU = 2 * mapHeight - 1, WU = 2 * mapWidth - 1;

	double**Umatr = new double*[HU];
	Umatr[0] = new double[HU*WU];
	for (int i = 1; i < HU; i++)
		Umatr[i] = Umatr[i - 1] + WU;

	for (int i = 0; i < HU; i++)
		for (int j = 0; j < WU; j++)
			Umatr[i][j] = 0;

	double dx = 0, dy = 0, dz1 = 0, dz2 = 0;
	if (cellType == RECT) {
		for (int i = 0; i < mapHeight; i++)
			for (int j = 0; j < mapWidth; j++)
			{
				if (j < mapWidth - 1)
				{
					for (int k = 0; k < variableNumber; k++)
						dx += pow(M[k][i][j] - M[k][i][j + 1], 2);
					Umatr[2 * i][2 * j + 1] = sqrt(dx);
				}
				if (i < mapHeight - 1)
				{
					for (int k = 0; k < variableNumber; k++)
						dy += pow(M[k][i][j] - M[k][i + 1][j], 2);
					Umatr[2 * i + 1][2 * j] = sqrt(dy);
				}
				if (i < (mapHeight - 1) && j < (mapWidth - 1))
				{
					for (int k = 0; k < variableNumber; k++)
					{
						dz1 += pow(M[k][i][j] - M[k][i + 1][j + 1], 2);
						dz2 += pow(M[k][i][j + 1] - M[k][i + 1][j], 2);
					}
					Umatr[2 * i + 1][2 * j + 1] = (sqrt(dz1) + sqrt(dz2)) / (2 * sqrt(2));
				}
				dx = dz1 = dz2 = dy = 0;
			}
	}
	else
	{
		for (int i = 0; i < mapHeight; i++)
			for (int j = 0; j < mapWidth; j++)
			{
				if (j < mapWidth - 1)
				{
					for (int k = 0; k < variableNumber; k++)
						dx += pow(M[k][i][j] - M[k][i][j + 1], 2);
					Umatr[2 * i][2 * j + 1] = sqrt(dx);
				}
				if (i < mapHeight - 1)
				{
					for (int k = 0; k < variableNumber; k++)
						dy += pow(M[k][i][j] - M[k][i + 1][j], 2);
					Umatr[2 * i + 1][2 * j] = sqrt(dy);
				}

				if ((i % 2 != 0) && i < (mapHeight - 1) && j < (mapWidth - 1))
				{
					for (int k = 0; k < variableNumber; k++)
						dz1 += pow(M[k][i][j] - M[k][i + 1][j + 1], 2);
					Umatr[2 * i + 1][2 * j + 1] = sqrt(dz1);
				}
				else if ((i % 2 == 0) && j > 0)
				{
					for (int k = 0; k < variableNumber; k++)
						dz1 += pow(M[k][i][j] - M[k][i + 1][j - 1], 2);
					Umatr[2 * i + 1][2 * j - 1] = sqrt(dz1);
				}
				dx = dz1 = dz2 = dy = 0;
			}
	}
	delete[] M[0][0];
	delete[]M[0];
	delete[] M;

	double mean = 0;
	int nX, nY, countNeighboor = 0;
	if (cellType == RECT) {
		for (int i = 0; i < HU; i += 2)
			for (int j = 0; j < WU; j += 2)
			{
				for (int k = -1; k < 2; k += 2)
				{
					nY = i + k;
					if (nY >= 0 && nY < HU)
					{
						mean += Umatr[nY][j];
						countNeighboor++;
					}
				}
				for (int k = -1; k < 2; k += 2)
				{
					nX = j + k;
					if (nX >= 0 && nX < WU)
					{
						mean += Umatr[i][nX];
						countNeighboor++;
					}
				}
				Umatr[i][j] = mean / countNeighboor;
				countNeighboor = 0;
				mean = 0;
			}
	}
	else
	{
		bool isOdd = false;
		for (int i = 0; i < HU; i += 2)
		{
			for (int j = 0; j < WU; j += 2)
			{
				for (int k = -1; k < 2; k += 2)
				{
					nY = i + k;
					if (nY >= 0 && nY < HU)
					{
						mean += Umatr[nY][j];
						countNeighboor++;
					}
				}
				for (int k = -1; k < 2; k += 2)
				{
					nX = j + k;
					if (nX >= 0 && nX < WU)
					{
						mean += Umatr[i][nX];
						countNeighboor++;
					}
				}

				if (isOdd)
				{
					nY = i - 1;
					nX = j + 1;
					if ((nX >= 0 && nX < WU) && (nY >= 0 && nY < HU))
					{
						mean += Umatr[nY][nX];
						countNeighboor++;
					}
					nY = i + 1;
					nX = j + 1;
					if ((nX >= 0 && nX < WU) && (nY >= 0 && nY < HU))
					{
						mean += Umatr[nY][nX];
						countNeighboor++;
					}
				}
				else
				{
					nY = i + 1;
					nX = j - 1;
					if ((nX >= 0 && nX < WU) && (nY >= 0 && nY < HU))
					{
						mean += Umatr[nY][nX];
						countNeighboor++;
					}
					nY = i - 1;
					nX = j - 1;
					if ((nX >= 0 && nX < WU) && (nY >= 0 && nY < HU))
					{
						mean += Umatr[nY][nX];
						countNeighboor++;
					}

				}

				Umatr[i][j] = mean / countNeighboor;
				countNeighboor = 0;
				mean = 0;
			}
			if (isOdd)isOdd = false;
			else isOdd = true;
		}
	}
	return Umatr;
}

void NewSOM::getNormalizeVector(const char* path)
{
	ifstream fin(path);
	double number = 0;
	for (int i = 0; i < variableNumber; i++)
	{
		fin >> number;
		normalizeVectorMAX[i] = number;
		normalizeVectorMIN[i] = number;
	}
	while (!fin.eof())
	{
		for (int i = 0; i < variableNumber; i++)
		{
			fin >> number;
			if (number > normalizeVectorMAX[i])
				normalizeVectorMAX[i] = number;
			if (number < normalizeVectorMIN[i])
				normalizeVectorMIN[i] = number;
		}
	}
	//bigSigma = 0;
	//for (int i = 0; i < variableNumber; i++)
	//	bigSigma += pow((normalizeVectorMAX[i] - normalizeVectorMIN[i])*sqrt(variableNumber), 2);
	//bigSigma = sqrt(bigSigma);
	string q = string(path);
	q = q.substr(0, q.length() - 4);
	if (clFile.is_open())clFile.close();
	clFile.open(q + "Clusters.txt", ios_base::out);

	if (qeHes.is_open())qeHes.close();
	qeHes.open(q + "qeHe.txt", ios_base::out);
	fin.close();
}

void NewSOM::saveVD()
{
	saveMap("123.cod");
	ofstream fv2("V.txt", ios_base::out);
	for (CLUSTERS::iterator it = clusters.begin(); it != clusters.end(); it++)
	{
		fv2 << it->first << ": ";
		for (set<int>::iterator ind = it->second.begin(); ind != it->second.end(); ind++)
		{
			fv2 << *ind << " ";
		}
		fv2<<endl;
	}
	fv2.close();

	ofstream fv("VV.txt", ios_base::out);
	for (int i = 0; i <2 * mapHeight - 1; i++)
	{
		for (int j = 0; j < 2 * mapWidth - 1; j++)
			if (i % 2 != 0 || j % 2 != 0)
			{
				if (V[i][j] <= 0)
					fv << setw(5) << right << -2<< ' ';
				else fv << setw(5) << right << -1 << ' ';
			}
			else
			{
				int neuronNum = (j / 2)*mapHeight + (i / 2);
				int a = 0;
				for (CLUSTERS::iterator it = clusters.begin(); it != clusters.end(); it++)
				{
					set<int>::iterator ind = it->second.find(neuronNum);
					if(ind!= it->second.end())
					{
						a++;
						fv << setw(5) << right << it->first << ' ';
						break;
					}
				}
				if(a==0)
					fv << setw(5) << right << -5 << ' ';
			}
			fv << endl;
	}
	fv.close();

	fv.open("VV3.txt", ios_base::out);
	for (int i = 0; i <2 * mapHeight - 1; i++)
	{
		for (int j = 0; j < 2 * mapWidth - 1; j++)
			if (i % 2 != 0 || j % 2 != 0)
			{
				if (V[i][j] <= 0)
					fv << setw(5) << right << '-' << ' ';
				else fv << setw(5) << right << '+' << ' ';
			}
			else
			{
				int neuronNum = (j / 2)*mapHeight + (i / 2);
				int a = 0;
				for (CLUSTERS::iterator it = clusters.begin(); it != clusters.end(); it++)
				{
					set<int>::iterator ind = it->second.find(neuronNum);
					if (ind != it->second.end())
					{
						a++;
						fv << setw(5) << right << it->first << ' ';
						break;
					}
				}
				if (a == 0)
					fv << setw(5) << right << '*' << ' ';
			}
		fv << endl;
	}
	fv.close();

	ofstream fv3("VV2.txt", ios_base::out);
	for (int i = 0; i <2 * mapHeight - 1; i++)
	{
		for (int j = 0; j < 2 * mapWidth - 1; j++)
			if (i % 2 != 0 || j % 2 != 0)
			{
				if (V[i][j] <= 0)
					fv3 << setw(5) << right << -1 << ' ';
				else fv3 << setw(5) << right << 1 << ' ';
			}
			else
			{
				int neuronNum = (j / 2)*mapHeight + (i / 2);
				fv3 << setw(5) << right << neuronNum << ' ';
			}
		fv3 << endl;
	}
	fv3.close();

	ofstream fd("VD.txt", ios_base::out);
	for (int i = 0; i < mapHeight; i++)
	{
		for (int j = 0; j < mapWidth; j++)
		{
			fd << D[mapHeight*j + i] << ' ';
		}
		fd << endl;
	}
	fd.close();

}

int** NewSOM::getAdjacencyMatrix()
{

	int** A = new int*[neuronNumber];
	int x, y, nX, nY,neuronN=0;

	A[0] = new int[neuronNumber * neighborsNumber];
	for (int i = 1; i < neuronNumber; i++)
		A[i] = A[i - 1] + neighborsNumber;

	for (int i = 0; i < neuronNumber; i++)
	{
		x = i / mapHeight;
		y = i%mapHeight;
		nY = y - 1;

		neuronN = x*mapHeight + nY;
		if (nY >= 0 && V[2 * y - 1][2 * x]>0)
		{
			A[i][0] = neuronN;
		}
		else { A[i][0] = -1; }

		nX = x + 1;

		neuronN = nX*mapHeight + y;
		if (nX < mapWidth && V[2 * y][2 * x + 1]>0)
		{
			A[i][1] = neuronN;
		}
		else { A[i][1] = -1; }
		nY = y + 1;
		neuronN = x*mapHeight + nY;
		if (nY < mapHeight && V[2 * y + 1][2 * x]>0)
		{
			A[i][2] = neuronN;
		}
		else { A[i][2] = -1; }
		nX = x - 1;
		neuronN = nX*mapHeight + y;
		if (nX >= 0 && V[2 * y][2 * x - 1]>0)
		{
			A[i][3] = neuronN;
		}
		else { A[i][3] = -1; }
		if (cellType == HEXA)
		{
			if (i % 2 != 0)
			{
				nY = y - 1;
				nX = x + 1;

				neuronN = nX*mapHeight + nY;
				if (nY >= 0 && nX < mapWidth && V[2 * y - 1][2 * x + 1]>0)
				{
					A[i][4] = neuronN ;
				}
				else { A[i][4] = -1; }

				nY = y + 1;
				nX = x + 1;
				neuronN = nX*mapHeight + nY;
				if (nY<mapHeight&&nX < mapWidth && V[2 * y + 1][2 * x + 1]>0)
				{
					A[i][5] = neuronN;
				}
				else { A[i][5] = -1; }
			}
			else
			{
				nY = y + 1;
				nX = x - 1;
				neuronN = nX*mapHeight + nY;
				if (nY<mapHeight && nX >0 && V[2 * y + 1][2 * x - 1]>0)
				{
					A[i][4] = neuronN;
				}
				else { A[i][4] = -1; }
				nY = y - 1;
				nX = x - 1;
				neuronN = nX*mapHeight + nY;
				if (nY >= 0 && nX >= 0 && V[2 * y - 1][2 * x - 1]>0)
				{
					A[i][5] = neuronN;
				}
				else { A[i][5] = -1; }
			}

		}
	}
	return A;
}

bool hasConnection(int*A,const int&len)
{
	for (int i = 0; i < len; i++)
	{
		if (A[i] != -1)return true;
	}
	return false;
}

void NewSOM::getClusters()
{
	int**AdM = getAdjacencyMatrix();

	int*group = new int[neuronNumber];
	setArray(group, neuronNumber, -2);
	int clustersN = 0;
	for (int i = 0; i < neuronNumber; i++)
	{
		if (hasConnection(AdM[i], neighborsNumber))
		{
			if (group[i] == -2)
			{
				group[i] = clustersN;
				findCluster(i, clustersN++, AdM,group);
			}
		}
		else
		{
			group[i] = -5;
		}
		if(group[i]!=-5)
			clusters[group[i]].insert(i);
	}

	delete[]AdM[0];
	delete[]AdM;
	delete[]group;
}
void NewSOM::findCluster(const int& ind,const int& cl,int**A,int*group)
{
	stack<int>neuronsToCheck;
	fillStack(neuronsToCheck, ind, A[ind],group);
	int newInd = 0;
	while(neuronsToCheck.size() != 0)
	{
		newInd =neuronsToCheck.top();
		neuronsToCheck.pop();
		if(group[newInd]==-2)
			fillStack(neuronsToCheck, newInd, A[newInd],group);
		group[newInd] = cl;
		
	}

}
void NewSOM::fillStack(stack<int>&neuronToCheck, const int&ind, int*A, int*group)
{
	//findNeighbors(ind);
	for (int i = 0; i < neighborsNumber; i++)
	{
		if (A[i] != -1)
		{
			if (group[A[i]] == -2)
			{
				neuronToCheck.push(A[i]);
			}
		}
	}
}

int NewSOM::getClustersNumber()
{
	return clusters.size();
}

void NewSOM::findSubClusters()
 {
	int max = getClustersNumber();
	vector<int> needToDelete;
	for (int i = 0; i != max; i++)
	{
		if (!subCluster(clusters[i], max,i))
		{
			needToDelete.push_back(i);
		}
	}
	deleteClusters(needToDelete);
	//int num = 0;
	//for (CLUSTERS::iterator it = clusters.begin(); it != clusters.end(); it++)
	//{
	//	if (num != it->first)
	//	{
	//		clusters[num] = clusters[it->first];
	//		needToDelete.push_back(it->first);
	//	}
	//	num++;
	//}
	//deleteClusters(needToDelete);
}

int findIndex(vector<int>&cluster,const int& val)
{
	for (int i = 0; i < cluster.size(); i++)
		if (val == cluster[i])
			return i;
	return -1;
}

bool NewSOM::subCluster(set<int>&cluster, const int& max,const int& I)
{
	vector<int> m;
	double**S;
	findMaxDensity(m,cluster);

	int maxNum = m.size();
	if (maxNum > 1)
	{
		CLUSTER newCluster;
		CLUSTERS cl;
		for (set<int>::iterator it = cluster.begin(); it != cluster.end(); it++)
		{
			newCluster[*it] = -2;
		}
		S = new double*[maxNum];
		S[0] = new double[maxNum*maxNum];
		for (int i = 1; i < maxNum; i++)
			S[i] = S[i - 1] + maxNum;

		//for (int i = 0; i < maxNum; i++)
		//{
		//	for (int j = 1; j < maxNum; j++)
		//	{
		//		if (D[m[i]] > (D[m[j]]))
		//		{
		//			m[j] = m[i] + m[j];
		//			m[i] = m[j] - m[i];
		//			m[j] = m[j] - m[i];
		//		}
		//	}
		//}



		for (int i = 0; i < maxNum; i++)
		{
			for (int j = 0; j < maxNum; j++)
			{
				S[i][j] = pow(((1 / D[m[i]]) + (1 / D[m[j]])), -1);
			}
			newCluster[m[i]] = m[i];
		}

		for (int j = 0; j < maxNum; j++)
		{
			makeSubCluster(m[j], newCluster);
		}

		//for (CLUSTER::iterator itI = newCluster.begin(); itI != newCluster.end(); itI++)
		//{
		//	for (CLUSTER::iterator itJ = newCluster.begin(); itJ != newCluster.end(); itJ++)
		//	{
		//		if (isNeuronsNeighbors(itI->first, itJ->first))
		//		{
		//			if (itI->second!=itJ->second)
		//			{
		//				int a = findIndex(m,itI->second),
		//				b = findIndex(m, itJ->second);

		//				if (!((D[itI->first] > S[a][b]) && (D[itJ->first] > S[a][b])))
		//				{
		//					if (D[itI->first] > D[itJ->first])
		//					{
		//						itJ->second = itI->second;
		//					}
		//					else
		//					{
		//						itI->second = itJ->second;
		//					}
		//				}
		//			}
		//		}
		//		
		//	}

		//}
		for (CLUSTER::iterator j = newCluster.begin(); j != newCluster.end(); j++)
		{
			if (j->second != -2)
			{
				cl[j->second].insert(j->first);
			}
		}

		for (int i = 0; i < maxNum; i++)
		{
			for (int j = i+1; j < maxNum; j++)
			{
				if (isClustersNeighbors(cl[m[i]], cl[m[j]])) 
				{
					if (isMerge(S, i, j, cl[m[i]], cl[m[j]]))
					{
						for (set<int>::iterator itI = cl[m[i]].begin(); itI != cl[m[i]].end(); itI++)
						{
							newCluster[*itI] = m[j];
							cl[m[j]].insert(*itI);
						}
						cl.erase(m[i]);
						i++;
						j=i;
					}
				}
			}
		}
		for (CLUSTER::iterator it = newCluster.begin(); it != newCluster.end(); it++)
		{
			if (it->second != -2)
				clusters[it->second + max].insert(it->first);
		}
		delete[]S[0];
		delete[]S;
		return false;
	}
	if (maxNum == 0)
	{
		return false;
	}
	return true;
}
void NewSOM::deleteClusters(vector<int>& needToDelete)
{
	for (int i = 0; i < needToDelete.size(); i++)
	{
		clusters.erase(needToDelete[i]);
	}
	needToDelete.clear();
}
void NewSOM::makeSubCluster(const int& num, CLUSTER&newCluster)
{
	stack<int> subCluster;
	titleNeigbors(num, newCluster, subCluster);
	set<int> usedNeurons;
	usedNeurons.insert(num);
	while (subCluster.size() != 0)
	{
		int a = subCluster.top();
		subCluster.pop();
		if (usedNeurons.find(a) == usedNeurons.end()) {
			titleNeigbors(a, newCluster, subCluster);
			newCluster[a] = num;
			usedNeurons.insert(a);
		}
	}

}

void NewSOM::titleNeigbors(const int&num, CLUSTER&cluster, stack<int>&subCluster)
{
	CLUSTER::iterator ind;
	findNeighbors(num);
	for (int i = 0; i < neighborsNumber; i++)
	{
		if (neighborhood[i] != -1)
		{
			ind = cluster.find(neighborhood[i]);
			if (ind != cluster.end())
			{
				if (D[num] > D[neighborhood[i]]) {
					subCluster.push(neighborhood[i]);
				}
			}
		}
	
	}
}

bool NewSOM::isMerge(double**S, const int& I, const int& J,set<int>& cl1, set<int>& cl2)
{
	if (cl1.size() == 0 || cl2.size() == 0)
		return false;
	for (set<int>::iterator itI = cl1.begin(); itI != cl1.end(); itI++)
	{
		for (set<int>::iterator itJ = cl2.begin(); itJ != cl2.end(); itJ++)
		{
			if (isNeuronsNeighbors(*itI, *itJ))
			{
				if (!((D[*itI] > S[I][J]) && (D[*itJ] > S[I][J])))
					return false;
			}
		}
	}
	return true;	
}


bool NewSOM::isClustersNeighbors(set<int>& cl1, set<int>& cl2)
{
	for (set<int>::iterator itI = cl1.begin(); itI != cl1.end(); itI++)
	{
		for (set<int>::iterator itJ = cl2.begin(); itJ != cl2.end(); itJ++)
		{
			if (isNeuronsNeighbors(*itI,*itJ))
			{
				return true;
			}
		}
	}
	return false;
}

bool NewSOM::isNeuronsNeighbors(const int& ind,const int& neigh)
{
	findNeighbors(ind);
	for (int i = 0; i < neighborsNumber; i++)
		if (neighborhood[i] == neigh)
		{
			return true;
		}
	return false;
}

void NewSOM::findMaxDensity(vector<int>& max, set<int>& cluster)
{
	for (set<int>::iterator it=cluster.cbegin();it!=cluster.end(); it++)
	{
		if (isDensityMax(*it))
			max.push_back(*it);
	}
}

bool NewSOM::isDensityMax(const int& ind)
{
	findNeighbors(ind);
	for (int i = 0; i < neighborsNumber; i++)
	{
		if (neighborhood[i] != -1)
		{
			if (D[neighborhood[i]] > D[ind])
			{
				return false;
			}
		}
	}
	return true;
}