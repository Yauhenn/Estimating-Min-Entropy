//LAst Tra-la-la

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <bitset>
#include <map>
#include <fstream>

using namespace std;

typedef unsigned char byte;

double mostCommonValue(int mas[], int countOfSymbols);

double collision(int mas[], int countOfSymbols, int unicSymbols);
double binarySearchForCollision(double lowerBound, int unicSymbols);
double calcRightSideForCollision(double p, int k);
double F(double q, int k);

double compression(int mas[], int countOfSymbols, int unicSymbols);
double binarySearchForCompression(double lowerBound, int unicSymbols, int countOfSymbols, int v);
double calcRightSideForCompression(double p, int k, int countOfSymbols, int v);
double G(double p, int countOfSymbols, int v);

double MultiMCW(int mas[], int countOfSymbols, int unicSymbols);
int mostCommon(int mas[], int a, int b, int unicSymbols);
double binarySearchForMultiMCW(int r, int N);
double calcRightSideForMultiMCW(double p, int r, int N);
int longestOnes(int correct[], int N);
double findX10(double p, int r);

double lagPrediction(int mas[], int countOfSymbols);

int main() {
	setlocale(LC_ALL, "Russian");
	
	int bytesInSource = 2097152; //1048576

	int bitsPerSymbol = 1;

	int bitsInSource = bytesInSource * 8;
	int *sourceInBits = new int[bitsInSource];
	
	FILE *f;
	f = fopen("aa_ab", "rb");

	int k = 0;
	unsigned int buffer;
	for (int i = 0; i<bytesInSource; i++)
	{
		fread(&buffer, 1, 1, f);
		byte myByte = static_cast<byte>(buffer);
		bitset<8> myBits(myByte);
		for (int i = 0, N = myBits.size(); i < N; ++i) {
			sourceInBits[k++] = myBits[N - i - 1];
		}
	}
	fclose(f);
	
	int countOfSymbols = bitsInSource / bitsPerSymbol;
	cout << "Count of symbols: " << countOfSymbols << endl;
	int *masOfSymbols = new int[countOfSymbols];

	if (bitsPerSymbol == 1) {
		for (int i = 0; i < countOfSymbols; i++)
			masOfSymbols[i] = sourceInBits[i];
	}
	else 
	{
		int k = 0;
		for (int i = 0; i < bitsInSource; i = i + bitsPerSymbol) {
			
			masOfSymbols[k] = 0;
			for (int j = 0; j < bitsPerSymbol; j++) {
				if (sourceInBits[i + j] == 1)
					masOfSymbols[k] += pow(2, bitsPerSymbol - 1 - j);
			}
			k++;
		}
	}

	delete[] sourceInBits;
	
	/*
	ofstream fout("cppstudio.txt");
	for (int i = 0; i < countOfSymbols; i++)
		fout << i << ": "<<masOfSymbols[i]<< endl;
	fout.close();
	*/
	int unicSymbols = pow(2, bitsPerSymbol);
	
	// Most Common Value
	cout << endl << "Most Common Value" << endl;
	double minEntropy1 = mostCommonValue(masOfSymbols, countOfSymbols);
	cout << "Entropy: " << minEntropy1 << endl;


	//Collision
	cout << endl << "Collision" << endl;
	double minEntropy2 = collision(masOfSymbols, countOfSymbols, unicSymbols);
	cout << "Entropy: " << minEntropy2 << endl;
	
	
	//Compression
	cout << endl << "Compression" << endl;
	double minEntropy3 = compression(masOfSymbols, countOfSymbols, unicSymbols);
	cout << "Entropy: " << minEntropy3 << endl;
	
	
	//MultiMCW
	cout << endl << "MultiMCW" << endl;
	double minEntropy4 = MultiMCW(masOfSymbols, countOfSymbols, unicSymbols);
	cout << "Entropy: " << minEntropy4 << endl;
	

	//Lag Prediction
	cout << endl << "Lag Prediction" << endl;
	double minEntropy5 = lagPrediction(masOfSymbols, countOfSymbols);
	cout << "Entropy: " << minEntropy5 << endl;

	delete[] masOfSymbols;

	system("pause");
	return 0;
}

double mostCommonValue(int mas[], int countOfSymbols) {
	map<int, int> dict;
	for (int i = 0; i < countOfSymbols; i++)
		dict[mas[i]]++;

	int maxCount = 0;
	map<int, int>::iterator it;
	for (it = dict.begin(); it != dict.end(); it++) {
		if (it->second >= maxCount)
			maxCount = it->second;
		
		cout << it->first << ": " << it->second << endl;
	}
	cout << "maxCount:" << maxCount << endl;

	double pMax = (double)maxCount / countOfSymbols;
	cout << "pMax: " << pMax << endl;

	double upperBound = pMax + 2.576 * sqrt(pMax * (1 - pMax) / countOfSymbols);
	cout << "upperBound: " << upperBound << endl;

	double pU = 1 < upperBound ? 1 : upperBound;
	cout << "pU: " << pU << endl;
	
	return -log10(pU)/log10(2);
}

double collision(int mas[], int countOfSymbols, int unicSymbols) {
	int v = 0;
	int index = 0;
	int *t = new int[countOfSymbols];

	for (int j = index; j < countOfSymbols; j++)
		for (int i = index; i < j; i++)
			if (mas[i] == mas[j]) {
				t[v] = j - index;
				v++;
				index = j + 1;
				break;
			}

	v--;
	cout << "v: " << v << endl;
	if (v < 1000) {
		cout << "Count of collisions is less than 1000, must take new data" << endl;
		return 0;
	}

	double meanX = 0.0;
	for (int i = 0; i < v; i++)
		meanX += t[i];
	meanX /= v;
	cout << "meanX: " << meanX << endl;

	double sigma = 0.0;
	for (int i = 0; i < v; i++)
		sigma += pow(t[i] - meanX, 2);
	delete[] t;
	sigma /= v;
	sigma = sqrt(sigma);
	cout << "sigma: " << sigma << endl;

	double lowerBound = meanX - 2.576 * sigma / sqrt(v);
	cout << "lowerBound: " << lowerBound << endl;

	double p = binarySearchForCollision(lowerBound, unicSymbols);
	cout << "p: " << p << endl;

	if (p > 0 && p < 1) {
		cout << "norm p" << endl;
		return -log10(p) / log10(2);
	}
	else{
		cout << "ne norm p" << endl;
		return log10(unicSymbols) / log10(2);
	}
}

double binarySearchForCollision(double lowerBound, int unicSymbols) {
	double minP = 1.0 / unicSymbols;
	double p_c = (1 - minP) / 2 + minP;
	double adj = 1 - minP;
	double rightSide = calcRightSideForCollision(p_c, unicSymbols);
	double maxTrueRightSide = calcRightSideForCollision(minP, unicSymbols);
	
	if (lowerBound > maxTrueRightSide) {
		cout << "something wrong!" << endl;
		return 0;
	}

	while (abs(lowerBound - rightSide) > 0.000001) {
		adj /= 2;
		if (lowerBound < rightSide) {
			p_c += adj;
			if (p_c == 1.0)
				p_c -= 0.000001;
		}
		else {
			p_c -= adj;
			if (p_c < minP)
				p_c = minP;
		}
		rightSide = calcRightSideForCollision(p_c, unicSymbols);
	}
	return p_c;
}

double calcRightSideForCollision(double p, int k) {
	double q = (1 - p) / (k - 1);
	double invQ = 1.0 / q;
	double invP = 1.0 / p;
	double invK = 1.0 / k;
	
	return p * pow(invQ, 2) * (1 + invK * (invP - invQ)) * F(q, k) - p * invQ * invK * (invP - invQ);
}

double F(double q, int k) {
	double z = 1 / q;
	double denominator = 1 + k / z;
	
	for (int i = 1; i < k; i++) {
		denominator = z + (-i) / denominator;
		denominator = 1 + (k - i) / denominator;
	}
	denominator = z + (-k) / denominator;

	return 1 / denominator;
}

double compression(int mas[], int countOfSymbols, int unicSymbols) {
	int d = 1000;
	int v = countOfSymbols - d;

	map<int, int> dict;
	for (int i = 0; i < unicSymbols; i++)
		dict[i] = 0;

	for (int i = 0; i < d; i++)
		dict[mas[i]] = i;

	int *D = new int[v];
	for (int i = d; i < countOfSymbols; i++) {
		if (dict[mas[i]] != 0) {
			D[i - d] = i - dict[mas[i]];
			dict[mas[i]] = i;
		}
		if (dict[mas[i]] == 0) {
			dict[mas[i]] = i;
			D[i - d] = i;
		}
	}

	int b = (int)(log10(unicSymbols) / log10(2) + 1);

	double meanX = 0.0;
	for (int i = 0; i < v; i++)
		meanX += log10(D[i]) / log10(2);
	meanX /= v;
	cout << "meanX: " << meanX << endl;

	double c = 0.7 - 0.8 / b + (4 + 32.0 / b) * pow(v, -3.0 / b) / 15;
	cout << "c: " << c << endl;

	double sigma = 0.0;
	for (int i = 0; i < v; i++)
		sigma += pow(log10(D[i]) / log10(2), 2);
	delete[] D;
	sigma /= v;
	sigma -= pow(meanX, 2);
	sigma = sqrt(sigma);
	sigma *= c;
	cout << "sigma: " << sigma << endl;

	double lowerBound = meanX - 2.576 * sigma / sqrt(v);
	cout << "lowerBound" << lowerBound << endl;

	double p = binarySearchForCompression(lowerBound, unicSymbols, countOfSymbols, v);
	cout << "p: " << p << endl;

	if (p > 0 && p < 1) {
		return -log10(p) / log10(2);
	}
	else
		return log10(unicSymbols) / log10(2);

}

double binarySearchForCompression(double lowerBound, int unicSymbols, int countOfSymbols, int v){
	double minP = 1.0 / unicSymbols;
	double p = (1 - minP) / 2 + minP;
	double adj = 1 - minP;
	double maxTrueRightSide = calcRightSideForCompression(minP, unicSymbols, countOfSymbols, v);
	if (lowerBound > maxTrueRightSide){
		cout << "something wrong!" << endl;
		return 0;
	}

	double rightSide = calcRightSideForCompression(p, unicSymbols, countOfSymbols, v);
	
	while (abs(lowerBound - rightSide) > 0.000001) {
		adj /= 2;
		if (lowerBound > rightSide)
			p -= adj;
		else
			p += adj;
		rightSide = calcRightSideForCompression(p, unicSymbols, countOfSymbols, v);
	}
	return p;
}

double calcRightSideForCompression(double p, int k, int countOfSymbols, int v) {
	double q = (1 - p) / (k - 1);
	return G(p, countOfSymbols, v) + (k - 1) * G(q, countOfSymbols, v);
}

double G(double p, int countOfSymbols, int v) {
	int d = countOfSymbols - v;
	double commonSum = 0.0;
	double innerSum = 0.0;
	for (int u = 1; u < d + 1; u++)
		innerSum += log10(u) / log10(2) * p * p *pow(1 - p, u - 1);
	commonSum = innerSum * v;

	double *another = new double[v];
	for (int t = d + 1; t < countOfSymbols + 1; t++)
		another[t - d - 1] = log10(t) / log10(2) * pow(1 - p, t - 1);

	double tempSum = 0.0;
	for (int i = 0; i < countOfSymbols - (d + 1); i++)
		tempSum += (countOfSymbols - i - (d + 1)) * another[i];
	commonSum += tempSum * p * p;

	double anotherSum = 0.0;
	for (int i = 0; i < v; i++)
		anotherSum += another[i];
	delete[] another;
	commonSum += p * anotherSum;
	
	return commonSum / v;

}

double MultiMCW(int mas[], int countOfSymbols, int unicSymbols) {
	int w[4] = {63, 255, 1023, 4095};
	int N = countOfSymbols - w[1];

	int *correct = new int[N];
	for (int i = 0; i < N; i++)
		correct[i] = 0;

	int scoreboard[4] = { 0, 0, 0, 0 };
	int frequent[4] = { -1, -1, -1, -1 };

	int winner = 0;

	for (int i = w[1]; i < countOfSymbols; i++) {
		if (i % 1000 == 0)
		cout << i << endl;
		for (int j = 0; j < 4; j++) {
			if (i > w[j])
				frequent[j] = mostCommon(mas, i - w[j], i - 1, unicSymbols);
			else
				frequent[j] = -1;
		}

		int prediction = frequent[winner];

		if (prediction == mas[i])
			correct[i - w[1]] = 1;

		for (int j = 0; j < 4; j++) {
			if (frequent[j] == mas[i]) {
				scoreboard[j] = scoreboard[j] + 1;
				if (scoreboard[j] >= scoreboard[winner])
					winner = j;
			}
		}

	}


	int C = 0;
	for (int i = 0; i < N; i++)
		if (correct[i] == 1)
			C++;

	double pGlobal = (float)C / N;
	cout << "pGlobal: "<<pGlobal << endl;
	
	double pGlobalUpperBound = pGlobal + 2.576 * sqrt(pGlobal * (1 - pGlobal) / (N - 1));
	cout << "pGlobalUpperBound: " << pGlobalUpperBound << endl;

	int r = longestOnes(correct, N) + 1;
	cout << "r: " << r << endl;
	delete[] correct;

	double pLocal = binarySearchForMultiMCW(r, N);
	cout << "pLocal: " << pLocal << endl;

	if (pGlobalUpperBound > pLocal)
		return -log10(pGlobalUpperBound) / log10(2);
	else
		return -log10(pLocal) / log10(2);
}

int mostCommon(int mas[], int a, int b, int unicSymbols) {
	map <int, int> dict;
	for (int i = a; i <= b; i++)
		dict[mas[i]]++;

	int maxCount = 0;
	int *list = new int[unicSymbols];
	int len = 0;

	map<int, int>::iterator it;
	for (it = dict.begin(); it != dict.end(); it++) {
		if (it->second == maxCount) {
			len++;
			list[len] = it->first;
		}
		if (it->second > maxCount) {
			len = 0;
			list[len] = it->first;
			maxCount = it->second;
		}
	}

	if (len == 0) 
		return list[len];
	else
	{
		for (int i = b; i >= a; i--)
			for (int j = 0; j < len; j++)
				if (mas[i] == list[j]) 
					return list[j];
	}
}

int longestOnes(int correct[], int N) {
	int run = 0;
	int longestRun = 0;
	for (int i = 0; i < N; i++) {
		if (correct[i] == 1)
			run++;
		else {
			if (longestRun < run)
				longestRun = run;
			run = 0;
		}
	}
	if (longestRun < run)
		longestRun = run;

	return longestRun;
}

double binarySearchForMultiMCW(int r, int N) {
	double p = 0.5;
	double adj = 0.5;
	double leftSide = 0.99;

	double rightSide = calcRightSideForMultiMCW(p, r, N);

	while (abs(rightSide - leftSide) > 0.000001) {
		adj /= 2;
		if (rightSide > leftSide)
			p += adj;
		else
			p -= adj;
		rightSide = calcRightSideForMultiMCW(p, r, N);
	}

	return p;
}

double calcRightSideForMultiMCW(double p, int r, int N) {
	double q = 1 - p;
	double x = findX10(p, r);

	return (1 - p * x) / ((r + 1 - r * x) * q) / pow(x, N + 1);

}

double findX10(double p, int r) {
	double q = 1 - p;
	double x = 1;

	for (int i = 0; i < 10; i++)
		x = 1 + q * pow(p, r) * pow(x, r + 1);

	return x;
}

double lagPrediction(int mas[], int countOfSymbols) {
	int D = 128;
	int N = countOfSymbols - 1;
	int *lag = new int[D];
	for (int i = 0; i < D; i++)
		lag[i] = -1;

	int *correct = new int[N];
	for (int i = 0; i < N; i++)
		correct[i] = 0;

	int *scoreboard = new int[D];
	for (int i = 0; i < D; i++)
		scoreboard[i] = 0;

	int winner = 0;
	
	for (int i = 1; i < countOfSymbols; i++) {
		for (int d = 0; d < D; d++) {
			if (d < i)
				lag[d] = mas[i - d - 1];
			else
				lag[d] = -1;
		}

		int prediction = lag[winner];

		if (prediction == mas[i])
			correct[i - 1] = 1;

		for (int d = 0; d < D; d++) {
			if (lag[d] == mas[i]) {
				scoreboard[d] = scoreboard[d] + 1;
				if (scoreboard[d] >= scoreboard[winner])
					winner = d;
			}
		}
	}

	int C = 0;
	for (int i = 0; i < N; i++)
		if (correct[i] == 1)
			C++;
	cout << "C: " << C << endl;

	double pGlobal = (float)C / N;
	cout << "pGlobal: " << pGlobal << endl;

	double pGlobalUpperBound = pGlobal + 2.576 * sqrt(pGlobal * (1 - pGlobal) / (N - 1));
	cout << "pGlobalUpperBound: " << pGlobalUpperBound << endl;

	int r = longestOnes(correct, N) + 1;
	cout << "r: " << r << endl;
	delete[] correct;

	double pLocal = binarySearchForMultiMCW(r, N);
	cout << "pLocal: " << pLocal << endl;

	if (pGlobalUpperBound > pLocal)
		return -log10(pGlobalUpperBound) / log10(2);
	else
		return -log10(pLocal) / log10(2);
}
