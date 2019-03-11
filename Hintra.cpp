#include<iostream>
#include<sstream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include<climits>
#include<cfloat>
#include<random>
#include<time.h>
#include<omp.h>

using namespace std;

double SumTwoLogs(double a, double b){
	double l2 = log(2);
	double shift = 0;
	double maxLog = log(DBL_MAX);
	double minLog = log(DBL_MIN);
	maxLog -= l2;
	if (std::isinf(a))
		a = DBL_MIN_EXP;
	if (std::isinf(b))
		b = DBL_MIN_EXP;
	if (a < minLog)
		shift = minLog - a;
	if (b + shift < minLog)
		shift = minLog - b;
	if (a + shift > maxLog)
		shift = maxLog - a;
	if (b + shift > maxLog)
		shift = maxLog - b;

	double s = exp(a + shift) + exp(b + shift);
	if (std::isinf(log(s)))
		cout << "inf" << endl;
	return log(s) - shift;
}


int SearchString(vector<string>& vector, string element){
	for (int i = 0; i < vector.size(); i++)
		if (vector[i] == element)
			return i;

	return -1;
}

bool BinarySearchString(vector<string>& vector, string element, int start, int end, int& pos){
	if (start > end){
		pos = start;
		return false;
	}
	int m = (start + end) / 2;
	if (vector[m] == element){
		pos = m;
		return true;
	}
	else if (vector[m] < element)
		return BinarySearchString(vector, element, m + 1, end, pos);
	else
		return BinarySearchString(vector, element, start, m - 1, pos);
}


double BinomialDensityLog(int a, int b, double p){
	if ((p == 0 && a == 0) || (p == 1 && b == 0))
		return 0;

	double choose = lgamma(b + 1) - lgamma(a + 1) - lgamma(b - a + 1);
	return choose + a * log(p) + (b - a) * log(1 - p);
}


int* NodeOrder(int* parentVector, int noNodes){
	int* ord = new int[noNodes];
	int size = 0;
	for (int i = 0; i < noNodes; i++)
		if (parentVector[i] == 0){
			ord[size] = i + 1;
			size++;
		}

	int cur = 0;
	while (size < noNodes){
		for (int i = 0; i < noNodes; i++)
			if (parentVector[i] == ord[cur])
				ord[size++] = i + 1;
		cur++;
	}

	return ord;
}


double* ToTheta(double* pis, int noNodes, int* parentVector, int* order){
	double* thetas = new double[noNodes + 1];

	for (int i = 0; i <= noNodes; i++)
		thetas[i] = pis[i];
	for (int i = noNodes - 1; i >= 0; i--)
		thetas[parentVector[order[i] - 1]] += thetas[order[i]];

	return thetas;
}


double MarginalGivenTree(int noNodes, int* parentVector, int* refs, int* vars, double precision, int* orderByLevels, int mutCount, double& maxProb){
	double* ys = new double[mutCount + 1];
	double parts = 1 / precision;
	int* placeHolders = new int[mutCount];
	for (int i = 0; i < mutCount; i++)
		placeHolders[i] = 0;

	double sum = DBL_MIN_EXP;
	maxProb = -1 * DBL_MAX;
	int counter = 0;
	while (true){
		ys[0] = (placeHolders[0] - 0) / parts;
		for (int i = 1; i < mutCount; i++)
			ys[i] = (placeHolders[i] - placeHolders[i - 1]) / parts;
		ys[mutCount] = (parts - placeHolders[mutCount - 1]) / parts;

		double *thetas = ToTheta(ys, mutCount, parentVector, orderByLevels);
		
		counter++;
		double prob = 0;
		int ind = 0;
		for (int i = 0; i < noNodes; i++){
			if (vars[i] == 0)
				continue;
			prob = prob + BinomialDensityLog(vars[i], vars[i] + refs[i], thetas[ind + 1] / 2);
			ind++;
		}

		sum = SumTwoLogs(sum, prob);
		if (prob > maxProb)
			maxProb = prob;
		//}
		delete[] thetas;

		int x;
		for (x = mutCount - 1; x >= 0; x--)
			if (placeHolders[x] < parts)
				break;

		//next
		if (x >= 0){
			int value = placeHolders[x] + 1;
			for (int i = x; i < mutCount; i++)
				placeHolders[i] = value;
		}
		else
			break;
	}


	delete[] ys;
	delete[] placeHolders;

	maxProb -= log(counter);
	return sum - log(counter);
}


double SearchForTreeMax(int* currentParentVector, int mutCount, int cayley, double* treePriors,
	int** allParentVects, double* maxProbs, int& maxL){
	int* bestParentVector = new int[mutCount];
	int ind;
	for (int i = 0; i < mutCount; i++)
		bestParentVector[i] = 0;

	double bestLikelihood = -1 * DBL_MAX, currentLikelihood, currentPrior;

	//int cayley = pow(mutCount + 1, mutCount - 1);
	for (int index = 0; index < cayley; index++)
	{
		for (int i = 0; i < mutCount; i++)
			currentParentVector[i] = allParentVects[index][i];

		currentLikelihood = treePriors[index] + maxProbs[index];// +TreePriorLog(currentParentVector, mutCount, 0.1);

		if (currentLikelihood > bestLikelihood){
			for (int i = 0; i < mutCount; i++)
				bestParentVector[i] = currentParentVector[i];
			bestLikelihood = currentLikelihood;
			maxL = index;
		}//if
	}//index

	for (int i = 0; i < mutCount; i++)
		currentParentVector[i] = bestParentVector[i];

	delete[] bestParentVector;

	return bestLikelihood;
}


double IntegrateOverTopologies(int cayley, double* treePriors,
	double* marginalProbs, double* maxProbs, double& maxLikelihood){

	double minDblExp = log(DBL_MIN);
	double maxDblExp = log(DBL_MAX);

	double sumLikelihood = 0, maxMarg = -1 * DBL_MAX;;

	//int cayley = pow(mutCount + 1, mutCount - 1);
	double* marglikelihoods = new double[cayley];
	maxLikelihood = -1 * DBL_MAX;
	double addedValue = 0, logCayley = log(cayley);
	for (int index = 0; index < cayley; index++)
	{
		marglikelihoods[index] = treePriors[index] + marginalProbs[index];
		double tempMaxLH = treePriors[index] + maxProbs[index];

		if (addedValue < minDblExp + logCayley - marglikelihoods[index])
			addedValue = minDblExp + logCayley - marglikelihoods[index];
		if (marglikelihoods[index] > maxMarg)
			maxMarg = marglikelihoods[index];
		if (tempMaxLH > maxLikelihood){
			maxLikelihood = tempMaxLH;
		}
		//delete[] orderByLevels;
	}//index

	if (maxMarg + addedValue > maxDblExp - logCayley)
		addedValue = maxDblExp - logCayley - maxMarg;

	for (int index = 0; index < cayley; index++)
		sumLikelihood += exp(marglikelihoods[index] + addedValue);

	delete[] marglikelihoods;

	return log(sumLikelihood) - addedValue;
}


struct CauseEffect{
	string cause;
	vector<string> effects;
	int causeInd;
	vector<int> effectInds;
};


void DFSParVec(int* mutSet, int size, int* parVector, int cur, vector<int> curCause, vector<CauseEffect>& ces){
	int ind;
	if (cur == 0)
		curCause.push_back(0);
	else{
		for (ind = 0; ind < curCause.size(); ind++)
			if (mutSet[cur - 1] < curCause[ind])
				break;
		curCause.insert(curCause.begin() + ind, mutSet[cur - 1]);
	}
	string strCause = "";
	for (int i = 0; i < curCause.size(); i++)
		strCause += to_string(curCause[i]) + ";";
	CauseEffect ace;
	ace.cause = strCause;
	bool hasEffect = false;
	for (int i = 0; i < size; i++)
		if (parVector[i] == cur){
			hasEffect = true;
			string strEffect = to_string(mutSet[i]);
			ace.effects.push_back(strEffect);
			DFSParVec(mutSet, size, parVector, i + 1, curCause, ces);
		}
	ces.push_back(ace);
}

void TreeToFactors(int noMuts, int size, bool* binaryProfile, int* parentVector, vector<CauseEffect>& ces){
	ces.clear();
	int* mutSet = new int[size];
	int index = 0;
	for (int i = 0; i < noMuts; i++)
		if (binaryProfile[i])
			mutSet[index++] = i + 1;

	// DFS
	vector<int> curCause;
	DFSParVec(mutSet, size, parentVector, 0, curCause, ces);

	delete[] mutSet;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){
	double precision = 0.1, epsilon = 1e-10, alpha=1;
	int numth = 4, nr = 0, nc = 0, EMIters = 100;
	string refFile, varFile, outtrees, /*outll, outcf,*/ outfacts, outpost, outmargf, outmaxf, pathToFile;
	bool shouldUpdate = false;

	if (argc > 1){
		string str(argv[1]);
		if (str == "-u"){
			shouldUpdate = true;
			if (argc != 8){
				cout << "Incorrect number of input arguments!" << endl;
				exit(1);
			}
		}
		else if (argc != 7){
			cout << "Incorrect number of input arguments!" << endl;
			exit(1);
		}

		int add = (shouldUpdate ? 1 : 0);
		pathToFile = argv[1 + add];
		refFile = pathToFile + ".Rcounts";
		cout << refFile << endl;
		varFile = pathToFile + ".Vcounts";
		cout << varFile << endl;
		outfacts = pathToFile + ".facts";
		cout << outfacts << endl;
		outpost = pathToFile + ".post";
		cout << outpost << endl;
		outtrees = pathToFile + ".trees";
		cout << outtrees << endl;
		/*outll = pathToFile + ".lltrace";
		cout << outll << endl;
		outcf = pathToFile + ".cftrace";
		cout << outcf << endl;*/
		outmargf = pathToFile + ".margs";
		cout << outmargf << endl;
		outmaxf = pathToFile + ".maxs";
		cout << outmaxf << endl;

		nr = atoi(argv[2 + add]);
		cout << nr << endl;
		nc = atoi(argv[3 + add]);
		cout << nc << endl;
		precision = atof(argv[4 + add]);
		cout << precision << endl;
		EMIters = atoi(argv[5 + add]);
		cout << EMIters << endl;
		numth = atoi(argv[6 + add]);
		cout << numth << endl;
	}
	else{
		cout << "Incorrect number of input arguments!" << endl;
		exit(1);
	}

	ofstream outputTrees, /*outputLL, outputCF,*/ outputFacts, outputPost;
	outputTrees.open(outtrees);
	//outputLL.open(outll);
	//outputCF.open(outcf);
	outputFacts.open(outfacts);
	outputPost.open(outpost);

	int** refs = new int*[nr];
	int** vars = new int*[nr];
	bool** binaryProfiles = new bool*[nr];
	int* mutCounts = new int[nr];

	//reading inputs
	cout << "Reading the read count data..." << endl;
	ifstream refin(refFile);
	ifstream varin(varFile);
	string rline, vline;
	int maxMutCount = 0;
	for (int row = 0; row < nr; row++){
		refs[row] = new int[nc];
		vars[row] = new int[nc];
		binaryProfiles[row] = new bool[nc];
		vars[row] = new int[nc];
		mutCounts[row] = 0;
		getline(refin, rline);
		getline(varin, vline);
		stringstream rs(rline);
		stringstream vs(vline);
		for (int col = 0; col < nc; col++){
			rs >> refs[row][col];
			vs >> vars[row][col];
			binaryProfiles[row][col] = (vars[row][col] > 0);
			mutCounts[row] += (vars[row][col] > 0 ? 1 : 0);
		}
		if (mutCounts[row] > maxMutCount)
			maxMutCount = mutCounts[row];
	}

	//Topologies
	int cayleys[] = { 1, 2, 9, 64, 625, 7776 };
	/*int* cayleys = new int[maxMutCount];
	for (int i = 0; i < maxMutCount; i++)
		cayleys[i] = pow(i + 2, i);*/
	int*** allParentVecs = new int**[maxMutCount];
	for (int k = 0; k < maxMutCount; k++){
		string strnc;
		ostringstream Convert;
		Convert << k + 1;
		strnc = Convert.str();
		string topfile = ".//PhylogenySet//ParentVectors_K-" + strnc + ".txt"; //******* SHOULD BE CUSTOMIZED TO YOUR DIRECTORY *******
		ifstream topin(topfile);
		allParentVecs[k] = new int*[cayleys[k]];
		for (int row = 0; row < cayleys[k]; row++){
			allParentVecs[k][row] = new int[k + 1];
			for (int col = 0; col <= k; col++)
				topin >> allParentVecs[k][row][col];
		}
	}

	// Combinations, Mapping, and Priors
	cout << "Identifying the mutation combinations in all patients..." << endl;
	vector<string> combinations;
	vector<double*> combPriors;
	vector<int> combSizes;
	int* combMapping = new int[nr];
	for (int i = 0; i < nr; i++){
		string code;
		int size = 0;
		for (int j = 0; j < nc; j++){
			code += (binaryProfiles[i][j] ? "1" : "0");
			if (binaryProfiles[i][j])
				size++;
		}
		bool found = false;
		for (int k = 0; k < combinations.size(); k++)
			if (combinations[k] == code){
				combMapping[i] = k;
				found = true;
				break;
			}
		if (!found){
			combinations.push_back(code);
			combMapping[i] = combinations.size() - 1;
			double* priors = new double[cayleys[size - 1]];
			for (int k = 0; k < cayleys[size - 1]; k++)
				priors[k] = -log(cayleys[size - 1]);
			combPriors.push_back(priors);
			combSizes.push_back(size);
		}
	}

	int* combCounts = new int[combinations.size()];
	for (int i = 0; i < combinations.size(); i++)
		combCounts[i] = 0;
	for (int i = 0; i < nr; i++)
		combCounts[combMapping[i]]++;

	cout << "Breaking all possible trees into factors..." << endl;
	vector<string> allCauses, allEffects;
	vector<CauseEffect>** factorsPerCombinationTopology = new vector<CauseEffect>*[combinations.size()];
	int* firstInstances = new int[combinations.size()];
	for (int c = 0; c < combinations.size(); c++){
		factorsPerCombinationTopology[c] = new vector<CauseEffect>[cayleys[combSizes[c] - 1]];
		int firstInstanc = 0;
		while (combMapping[firstInstanc++] != c);
		firstInstanc--;
		firstInstances[c] = firstInstanc;
		for (int p = 0; p < cayleys[combSizes[c] - 1]; p++){
			TreeToFactors(nc, mutCounts[firstInstanc], binaryProfiles[firstInstanc], allParentVecs[mutCounts[firstInstanc] - 1][p], factorsPerCombinationTopology[c][p]);
			for (int f = 0; f < factorsPerCombinationTopology[c][p].size(); f++){
				int ind = SearchString(allCauses, factorsPerCombinationTopology[c][p][f].cause);
				if (ind == -1){
					allCauses.push_back(factorsPerCombinationTopology[c][p][f].cause);
					factorsPerCombinationTopology[c][p][f].causeInd = allCauses.size() - 1;
				}
				else {
					factorsPerCombinationTopology[c][p][f].causeInd = ind;
				}

				for (int e = 0; e < factorsPerCombinationTopology[c][p][f].effects.size(); e++){
					ind = SearchString(allEffects, factorsPerCombinationTopology[c][p][f].effects[e]);
					if (ind == -1){
						allEffects.push_back(factorsPerCombinationTopology[c][p][f].effects[e]);
						factorsPerCombinationTopology[c][p][f].effectInds.push_back(allEffects.size() - 1);
					}
					else {
						factorsPerCombinationTopology[c][p][f].effectInds.push_back(ind);
					}
				}
			}
		}
	}

	// Computing Individual Marginals
	cout << "Computing the marginal probabilities..." << endl;
	double** marginalProbs = new double*[nr];
	double** maxProbs = new double*[nr];
	ifstream inmargs(outmargf);
	ifstream inmaxs(outmaxf);
	if (shouldUpdate || !inmargs.good() || !inmaxs.good()){
		ofstream outmargs(outmargf);
		ofstream outmaxs(outmaxf);
		for (int r = 0; r < nr; r++){
			marginalProbs[r] = new double[cayleys[mutCounts[r] - 1]];
			maxProbs[r] = new double[cayleys[mutCounts[r] - 1]];
#pragma omp parallel for num_threads(numth)
			for (int index = 0; index < cayleys[mutCounts[r] - 1]; index++)
			{
				int* currentParentVector = new int[mutCounts[r]];
				int* orderByLevels;
				cout << "sample " << r << " tree " << index << endl;

				for (int i = 0; i < mutCounts[r]; i++)
					currentParentVector[i] = allParentVecs[mutCounts[r] - 1][index][i];
				orderByLevels = NodeOrder(currentParentVector, mutCounts[r]);
				marginalProbs[r][index] = MarginalGivenTree(nc, currentParentVector, refs[r], vars[r], precision, orderByLevels, mutCounts[r], maxProbs[r][index]);
				delete[] orderByLevels;
				delete[] currentParentVector;
			}//index

			for (int index = 0; index < cayleys[mutCounts[r] - 1]; index++)
			{
				outmargs << marginalProbs[r][index] << "\t";
				outmaxs << maxProbs[r][index] << "\t";
			}
			outmargs << endl;
			outmaxs << endl;
		}//r
	}
	else{
		for (int r = 0; r < nr; r++){
			marginalProbs[r] = new double[cayleys[mutCounts[r] - 1]];
			maxProbs[r] = new double[cayleys[mutCounts[r] - 1]];
			//#pragma omp parallel for num_threads(numth)
			for (int index = 0; index < cayleys[mutCounts[r] - 1]; index++)
			{
				inmargs >> marginalProbs[r][index];
				inmaxs >> maxProbs[r][index];
			}//index
		}
	}

	int** currentParentVectors = new int*[nr];
	for (int i = 0; i < nr; i++)
		currentParentVectors[i] = new int[mutCounts[i]];

	double currentMargLogProb = 0, currentMaxLogProb = 0, treeMargL = 0, treeMaxL = 0;
	double* currentMargLogProbs = new double[nr];
	double* currentMaxLogProbs = new double[nr];
	for (int i = 0; i < nr; i++){
		currentMargLogProbs[i] = IntegrateOverTopologies(cayleys[mutCounts[i]-1], combPriors[combMapping[i]], marginalProbs[i], maxProbs[i], currentMaxLogProbs[i]);
		currentMargLogProb += currentMargLogProbs[i];
		currentMaxLogProb += currentMaxLogProbs[i];

		int itemp;
		SearchForTreeMax(currentParentVectors[i], mutCounts[i], cayleys[mutCounts[i] - 1], combPriors[combMapping[i]], allParentVecs[mutCounts[i] - 1], marginalProbs[i], itemp);
		treeMaxL += maxProbs[i][itemp];
		treeMargL += marginalProbs[i][itemp];
	}

	double prevMargLogProb = currentMargLogProb;

	// Global Factor Priors
	double** factorPriors = new double*[allCauses.size()];
	//double** factorPriorsBar = new double*[allCauses.size()];
	double** proposedFactorPriors = new double*[allCauses.size()];
	//double** proposedFactorPriorsBar = new double*[allCauses.size()];
	double* causeCounts = new double[allCauses.size()];
	double* discCauseCounts = new double[allCauses.size()];
	double* proposedCauseCounts = new double[allCauses.size()];
	double* proposedDiscCauseCounts = new double[allCauses.size()];
	for (int i = 0; i < allCauses.size(); i++){
		factorPriors[i] = new double[allEffects.size()];
		//factorPriorsBar[i] = new double[allEffects.size()];
		proposedFactorPriors[i] = new double[allEffects.size()];
		//proposedFactorPriorsBar[i] = new double[allEffects.size()];
	}

	for (int c = 0; c < allCauses.size(); c++){
		causeCounts[c] = log(2 * epsilon);
		discCauseCounts[c] = log(epsilon);
		for (int e = 0; e < allEffects.size(); e++){
			factorPriors[c][e] = log(epsilon);
		}
	}

	cout << "Learning the tree parameters..." << endl;
	for (int trial = 0; trial < EMIters; trial++){
		// Reporting the current status
		cout << trial << " (" << currentMargLogProb << ")" << endl;
		/*outputLL << currentMargLogProb << "\t" << currentMaxLogProb << "\t" << treeMargL << "\t" << treeMaxL << endl;

		for (int c = 0; c < combinations.size(); c++){
			for (int index = 0; index < cayleys[combSizes[c] - 1]; index++)
				outputCF << combPriors[c][index] << "\t";
			outputCF << endl;
		}
		outputCF << endl;*/

		//// Updating Counts

		// Resetting the Factor Priors
		for (int c = 0; c < allCauses.size(); c++){
			proposedCauseCounts[c] = log(2 * epsilon);
			proposedDiscCauseCounts[c] = log(epsilon);
			for (int e = 0; e < allEffects.size(); e++){
				proposedFactorPriors[c][e] = log(epsilon);
				//proposedFactorPriorsBar[c][e] = log(2 * epsilon);
			}
		}

		// #pragma omp parallel for num_threads(numth)
		for (int i = 0; i < nr; i++){
			for (int t = 0; t < cayleys[mutCounts[i] - 1]; t++){
				vector<string> strUniqueCauses;
				vector<int> intUniqueCauses;
				double coef = combPriors[combMapping[i]][t] + marginalProbs[i][t] - currentMargLogProbs[i]; //confidence

				for (int f = 0; f < factorsPerCombinationTopology[combMapping[i]][t].size(); f++){
					proposedCauseCounts[factorsPerCombinationTopology[combMapping[i]][t][f].causeInd] =
						SumTwoLogs(coef, proposedCauseCounts[factorsPerCombinationTopology[combMapping[i]][t][f].causeInd]);
					if (factorsPerCombinationTopology[combMapping[i]][t][f].effects.size()==0)
						proposedDiscCauseCounts[factorsPerCombinationTopology[combMapping[i]][t][f].causeInd] =
						SumTwoLogs(coef, proposedDiscCauseCounts[factorsPerCombinationTopology[combMapping[i]][t][f].causeInd]);
					/*for (int j = 0; j < allEffects.size(); j++){
						int effind = stoi(allEffects[j]) - 1;
						if (binaryProfiles[firstInstances[combMapping[i]]][effind])
							proposedFactorPriorsBar[factorsPerCombinationTopology[combMapping[i]][t][f].causeInd][j] =
							SumTwoLogs(coef, proposedFactorPriorsBar[factorsPerCombinationTopology[combMapping[i]][t][f].causeInd][j]);
					}*/
					for (int e = 0; e < factorsPerCombinationTopology[combMapping[i]][t][f].effects.size(); e++){
						proposedFactorPriors[factorsPerCombinationTopology[combMapping[i]][t][f].causeInd][factorsPerCombinationTopology[combMapping[i]][t][f].effectInds[e]] =
							SumTwoLogs(coef, proposedFactorPriors[factorsPerCombinationTopology[combMapping[i]][t][f].causeInd][factorsPerCombinationTopology[combMapping[i]][t][f].effectInds[e]]);
					}
				}
			}
		} // i < nr

		//// Normalizing the Factor Counts
		//for (int c = 0; c < allCauses.size(); c++){
		//	for (int e = 0; e < allEffects.size(); e++)
		//		proposedFactorPriors[c][e] -= proposedFactorPriorsBar[c][e];
		//}

		// Updating the Combination Priors
#pragma omp parallel for num_threads(numth)
		for (int c = 0; c < combinations.size(); c++){
			double priorSum = log(DBL_MIN);
			for (int t = 0; t < cayleys[combSizes[c] - 1]; t++){
				double prior = 0;
				for (int f = 0; f < factorsPerCombinationTopology[c][t].size(); f++){
					/*for (int j = 0; j < allEffects.size(); j++){
						int effind = stoi(allEffects[j]) - 1;
						if (binaryProfiles[firstInstances[c]][effind])
							prior += log(1 - exp(proposedFactorPriors[factorsPerCombinationTopology[c][t][f].causeInd][j]));
					}*/
					for (int e = 0; e < factorsPerCombinationTopology[c][t][f].effects.size(); e++)
						prior += proposedFactorPriors[factorsPerCombinationTopology[c][t][f].causeInd][factorsPerCombinationTopology[c][t][f].effectInds[e]];// -
						//log(1 - exp(proposedFactorPriors[factorsPerCombinationTopology[c][t][f].causeInd][factorsPerCombinationTopology[c][t][f].effectInds[e]]));
				}

				combPriors[c][t] = prior*alpha;
				priorSum = SumTwoLogs(prior*alpha, priorSum);
			}// index

			for (int t = 0; t < cayleys[combSizes[c] - 1]; t++)
				combPriors[c][t] = combPriors[c][t] - priorSum;
		}// c

		currentMargLogProb = currentMaxLogProb = treeMaxL = treeMargL = 0;
		for (int i = 0; i < nr; i++){
			currentMargLogProbs[i] = IntegrateOverTopologies(cayleys[mutCounts[i] - 1], combPriors[combMapping[i]], marginalProbs[i], maxProbs[i], currentMaxLogProbs[i]);
			currentMargLogProb += currentMargLogProbs[i];
			currentMaxLogProb += currentMaxLogProbs[i];

			int itemp;
			SearchForTreeMax(currentParentVectors[i], mutCounts[i], cayleys[mutCounts[i] - 1], combPriors[combMapping[i]], allParentVecs[mutCounts[i] - 1], marginalProbs[i], itemp);
			treeMaxL += maxProbs[i][itemp];
			treeMargL += marginalProbs[i][itemp];
		}

		// Not improved?
		if (currentMargLogProb < prevMargLogProb){
		//if (false){
			// Updating the Combination Priors
#pragma omp parallel for num_threads(numth)
			for (int c = 0; c < combinations.size(); c++){
				double priorSum = log(DBL_MIN);
				for (int t = 0; t < cayleys[combSizes[c] - 1]; t++){
					double prior = 0;
					for (int f = 0; f < factorsPerCombinationTopology[c][t].size(); f++){
						/*for (int j = 0; j < allEffects.size(); j++){
							int effind = stoi(allEffects[j]) - 1;
							if (binaryProfiles[firstInstances[c]][effind])
								prior += log(1 - exp(factorPriors[factorsPerCombinationTopology[c][t][f].causeInd][j]));
						}*/
						for (int e = 0; e < factorsPerCombinationTopology[c][t][f].effects.size(); e++)
							prior += factorPriors[factorsPerCombinationTopology[c][t][f].causeInd][factorsPerCombinationTopology[c][t][f].effectInds[e]];// -
							//log(1 - exp(factorPriors[factorsPerCombinationTopology[c][t][f].causeInd][factorsPerCombinationTopology[c][t][f].effectInds[e]]));
					}

					combPriors[c][t] = prior*alpha;
					priorSum = SumTwoLogs(prior*alpha, priorSum);
				}// index

				for (int t = 0; t < cayleys[combSizes[c] - 1]; t++)
					combPriors[c][t] = combPriors[c][t] - priorSum;
			}// c

			for (int c = 0; c < allCauses.size(); c++){
				causeCounts[c] = proposedCauseCounts[c];
				discCauseCounts[c] = proposedDiscCauseCounts[c];
				for (int e = 0; e < allEffects.size(); e++){
					factorPriors[c][e] = proposedFactorPriors[c][e];
					//factorPriorsBar[c][e] = proposedFactorPriorsBar[c][e];
				}
			}
			break;
		}// if
		
		prevMargLogProb = currentMargLogProb;
		for (int c = 0; c < allCauses.size(); c++){
			causeCounts[c] = proposedCauseCounts[c];
			discCauseCounts[c] = proposedDiscCauseCounts[c];
			for (int e = 0; e < allEffects.size(); e++){
				factorPriors[c][e] = proposedFactorPriors[c][e];
				//factorPriorsBar[c][e] = proposedFactorPriorsBar[c][e];
			}
		}
	}// EM Iterations


	//// Outputs

	// Output Factor Parameters
	outputFacts << "Ancestors\tSupport\tSupport (Discont'd)\t";
	//outputTotFacts << "Ancestors\tSupport\t";
	for (int e = 0; e < allEffects.size(); e++){
		outputFacts << allEffects[e] << "\t";
		//outputTotFacts << allEffects[e] << "\t";
	}
	outputFacts << endl;
	//outputTotFacts << endl;
	for (int ca = 0; ca < allCauses.size(); ca++){
		outputFacts << allCauses[ca] << "\t" << exp(causeCounts[ca]) << "\t" << exp(discCauseCounts[ca]) << "\t";
		//outputTotFacts << allCauses[ca] << "\t" << exp(causeCounts[ca]) << "\t";
		for (int e = 0; e < allEffects.size(); e++){
			outputFacts << exp(factorPriors[ca][e]) << "\t";
			//outputTotFacts << exp(factorPriors[ca][e] + factorPriorsBar[ca][e]) << "\t";
		}
		outputFacts << endl;
		//outputTotFacts << endl;
	}

	// Compute the trees & count factors
	cout << "Computing the maximum a posteriori trees..." << endl;

	// Resetting the Factor Priors (not in log scale)
	for (int c = 0; c < allCauses.size(); c++){
		causeCounts[c] = 0;
		for (int e = 0; e < allEffects.size(); e++){
			factorPriors[c][e] = 0;
		}
	}

	int* parentVector = new int[nc];
	int** finalParentVectors = new int*[nr];
	for (int k = 0; k < nr; k++){
		int itemp;
		SearchForTreeMax(currentParentVectors[k], mutCounts[k], cayleys[mutCounts[k] - 1], combPriors[combMapping[k]], allParentVecs[mutCounts[k] - 1], marginalProbs[k], itemp);

		for (int f = 0; f < factorsPerCombinationTopology[combMapping[k]][itemp].size(); f++){
			causeCounts[factorsPerCombinationTopology[combMapping[k]][itemp][f].causeInd]++;
			for (int e = 0; e < factorsPerCombinationTopology[combMapping[k]][itemp][f].effects.size(); e++)
				factorPriors[factorsPerCombinationTopology[combMapping[k]][itemp][f].causeInd][factorsPerCombinationTopology[combMapping[k]][itemp][f].effectInds[e]]++;
		}


		int ind = 0;
		for (int i = 0; i < nc; i++)
			parentVector[i] = (binaryProfiles[k][i] ? currentParentVectors[k][ind++] : -1);
		for (int i = 0; i < nc; i++)
			if (!binaryProfiles[k][i])
				for (int j = 0; j < nc; j++)
					if (parentVector[j]>i)
						parentVector[j]++;
		finalParentVectors[k] = new int[nc];
		for (int i = 0; i < nc; i++)
			finalParentVectors[k][i] = parentVector[i];
	}

	// Output the trees
	for (int i = 0; i < nr; i++){
		for (int j = 0; j < nc; j++)
			outputTrees << finalParentVectors[i][j] << "\t";
		outputTrees << endl;
	}


	// Output Post Factor Parameters
	outputPost << "Ancestors\tSupport\t";
	for (int e = 0; e < allEffects.size(); e++){
		outputPost << allEffects[e] << "\t";
	}
	outputPost << endl;
	for (int ca = 0; ca < allCauses.size(); ca++){
		outputPost << allCauses[ca] << "\t" << causeCounts[ca] << "\t";
		for (int e = 0; e < allEffects.size(); e++){
			outputPost << factorPriors[ca][e] << "\t";
		}
		outputPost << endl;
	}
}
