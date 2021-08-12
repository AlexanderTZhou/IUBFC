/***
AUTHOR: ALEXANDER ZHOU
INSTITUTION: HONG KONG UNIVERSITY OF SCIENCE AND TECHNOLOGY
E-MAIL: atzhou@cse.ust.hk

PLEASE E-MAIL ANY QUERIES, BUT DEPENDING ON WHEN IT'S SENT THERE'S
A NON-ZERO (VERY LIKELY) CHANCE I'VE FORGOTTEN EVERYTHING ABOUT THIS CODE
***/

#define _CRT_SECURE_NO_WARNINGS

#define ll long long
#define ull unsigned long long
#define ld long double

#include<iostream>
#include<fstream>
#include<string>
#include<cstring>
#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<unordered_map>
#include<utility>
#include<functional>
#include<cstdio>
#include<ctime>
#include<chrono>
#include<algorithm>
#include<cmath>

using namespace std;


/***
* The node class.
* id = ID of the Node
* neighbours = all the neighbours of the node and their respective probability
*/
class Node {
	public:
		int id;
		vector<pair<int, double> > neighbours;
	
		//Inserts a neighbour (ID, Prob)
		void insertNeighbour(int n, double d) {
			neighbours.push_back(make_pair(n, d));
		}
};

/***
* The edge class (used predominately for pruning/sampling)
* n1, n2: ID of Node 1 and Node 2 respectively, n1 < n2
* prob: The existential probability of the edge
*/
class Edge {
	public:
		int n1;
		int n2;
		double prob;
};

/***
* A struct that may be used to compare the first element in
* a pair<int, double> to an int. Useful for Binary Search
*/
struct PairComp { 
	bool operator()(const pair<int, double>& v1, const pair<int, double>& v2) {
		return (v1.first < v2.first);
	}
    bool operator()(const pair<int, double>& value, const int& key) { 
        return (value.first < key); 
    } 
    bool operator()(const int& key, const pair<int, double>& value) { 
        return (key < value.first); 
    } 
}; 

vector<int> vID; //All IDs
unordered_map<int, Node> nodeMap; //ID and the actual Node
vector<Edge> allEdges; //Edges
double t, maxT; //The Threshold Probability and maximum threshold among edges

/***
* Initialises ID File (tsv format)
*/
void initialiseID(char* filename) {
	ifstream nameStream(filename);
	if (!nameStream.is_open()) {
		printf("Invalid File (vertices)\n");
		exit(2);
	}
	vector<int> vecName;
	char line[20];
	char* tempV;
	int id;

	while (nameStream.good()) {
		nameStream.getline(line, 20);
		if (nameStream.eof()) {
			nameStream.close();
			break;
		}
		tempV = strtok(line, "\t");
		id = stoi(tempV);
		Node n;
		n.id = id;
		vID.push_back(id);
		nodeMap.insert(make_pair(id, n));
	}
//	cout << "Nodes: " << vID.size() << "\n";
}

/***
* Initialise all edges, removing those that break bipartite format (tsv format)
*/

void initialiseEdge(char* filename) {
	ifstream nameStream(filename);
	if (!nameStream.is_open()) {
		printf("Edges File DNE\n");
		exit(3);
	}
	char line[40];
	char* tempV;
	int v1, v2;
	double p;
	int edgeCount = 0;
	maxT = 0;
	unordered_map<int, Node>::iterator nIt1, nIt2;
	while (nameStream.good()) {
		nameStream.getline(line, 40);
		edgeCount++;
		if (nameStream.eof()) {
			nameStream.close();
			break;
		}
		
		tempV = strtok(line, "\t");
		v1 = stoi(tempV);
		tempV = strtok(NULL, "\t");
		v2 = stoi(tempV);
		tempV = strtok(NULL, "\t");
		p = stod(tempV);
		
		nIt1 = nodeMap.find(v1);
		if (nIt1 == nodeMap.end()) {
			printf("Incorrect Edge File\n");
			exit(4);
		}
		
		nIt2 = nodeMap.find(v2);
		if (nIt2 == nodeMap.end()) {
			printf("Incorrect Edge File\n");
			exit(4);
		}
		
		nIt1->second.neighbours.push_back(make_pair(v2, p));
		nIt2->second.neighbours.push_back(make_pair(v1, p));
		
		Edge e;
		e.n1 = v1;
		e.n2 = v2;
		e.prob = p;
		allEdges.push_back(e);
		
		if (p > maxT) {
			maxT = p;
		}
	}
	
	//Turn on if EDFEFILE format is unverified/unsorted
	PairComp pc;
	for (auto &a : nodeMap) {
		sort(a.second.neighbours.begin(), a.second.neighbours.end(), pc);
	}
//	cout << "Edges: " << edgeCount << "\n";
}

/***
* Initialise all edges, removing those that break bipartite format (tsv format)
* Also prunes edges less than t
*/

void initialiseEdgePrune(char* filename) {
	ifstream nameStream(filename);
	if (!nameStream.is_open()) {
		printf("Edges File DNE\n");
		exit(3);
	}
	char line[30];
	char* tempV;
	int v1, v2;
	double p;
	int edgeCount = 0;
	int totEdge = 0;
	maxT = 0;
	unordered_map<int, Node>::iterator nIt1, nIt2;
	while (nameStream.good()) {
		
		totEdge++;
		nameStream.getline(line, 30);
		if (nameStream.eof()) {
			nameStream.close();
			break;
		}
		
		tempV = strtok(line, "\t");
		v1 = stoi(tempV);
		tempV = strtok(NULL, "\t");
		v2 = stoi(tempV);
		tempV = strtok(NULL, "\t");
		p = stod(tempV);
		if (p < t) {
			continue;
		}
		nIt1 = nodeMap.find(v1);
		if (nIt1 == nodeMap.end()) {
			printf("Incorrect Edge File\n");
			exit(4);
		}
		
		nIt2 = nodeMap.find(v2);
		if (nIt2 == nodeMap.end()) {
			printf("Incorrect Edge File\n");
			exit(4);
		}
		
		edgeCount++;
		
		nIt1->second.neighbours.push_back(make_pair(v2, p));
		nIt2->second.neighbours.push_back(make_pair(v1, p));
		
		Edge e;
		e.n1 = v1;
		e.n2 = v2;
		e.prob = p;
		allEdges.push_back(e);
		
		if (p > maxT) {
			maxT = p;
		}
	}
	
	//Turn on if EDFEFILE format is unverified/unsorted
	PairComp pc;
	for (auto &a : nodeMap) {
		sort(a.second.neighbours.begin(), a.second.neighbours.end(), pc);
	}
//	cout << "Edges: " << totEdge << "\n";
//	cout << "Deleted Edges: " << totEdge - edgeCount << "\n";
}

void writeToFile(string s, char* filename) {
	ofstream writeFile;
	writeFile.open(filename, ios::app);
	if (writeFile.is_open()) {
		writeFile << s << "\n";
		writeFile.flush();
	} 
	writeFile.close();
}


/***
* Operator for Vertex Priority of Two Nodes
*/
bool vpComp(Node n1, Node n2) {
	if (n1.neighbours.size() == n2.neighbours.size()) {
		return n1.id > n2.id;
	}
	return n1.neighbours.size() > n2.neighbours.size();
}

/***
* Reverse Operator for Vertex Priority of Two Nodes
* Used for the purposes of sorting
*/
bool vpRevComp(Node n1, Node n2) {
	if (n1.neighbours.size() == n2.neighbours.size()) {
		return n1.id < n2.id;
	}
	return n1.neighbours.size() < n2.neighbours.size();
}

/***
* Sorts the neighbours by (inc) Vertex Priority
*/
void vpNeighbourSort(Node &n) {	
	sort(n.neighbours.begin(), n.neighbours.end(), [](const pair<int, double> &p1, const pair<int, double> &p2) -> bool {
		return vpRevComp(nodeMap.find(p1.first)->second,nodeMap.find(p2.first)->second);
	});
}

/***
* Initialise the Vertex Priority (e.g. degreeness) of each node
*/
void initialiseVP() {
	for (unordered_map<int, Node>::iterator nIt = nodeMap.begin(); nIt != nodeMap.end(); nIt++) {
		vpNeighbourSort(nIt->second);
	}
}

/***
* Sorts the neighbours of a node by (dec) Existential Probability
*/
void epNeighbourSort(Node &n) {
	sort(n.neighbours.begin(), n.neighbours.end(), [](const pair<int, double> &p1, const pair<int, double> &p2) -> bool {
		return p1.second > p2.second;
	});
}

/***
* Initialise the Existential Probability sorted order of each node
*/
void initialiseEP() {
	for (unordered_map<int, Node>::iterator nIt = nodeMap.begin(); nIt != nodeMap.end(); nIt++) {
		epNeighbourSort(nIt->second);
	}
}

/***
* Finds and Calculates the butterfly existential probability
* for the baseline method
*/
double baseListProbCalc(int u, int w, vector<pair<int, double> > &edgeList){
	bool uTrig = false;
	bool wTrig = false;
	double uEProb, wEProb;
	for (pair<int, double> tp : edgeList) {
		if (tp.first == u) {
			uEProb = tp.second;
			uTrig = true;
		} else if (tp.first == w) {
			wEProb = tp.second;
			wTrig = true;
		}
		if (uTrig && wTrig) {
			return uEProb * wEProb;
		}
	}
	return 0.0;
}

/***
* Baseline method of counting the number of uncertain butterflies from a list
*/
long baseListCount(int u, int w, vector<int> &wedges) {
	long count = 0;
	int n = wedges.size();
	double wedgeProb1, wedgeProb2;
	unordered_map<int, Node>::iterator nIt;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			nIt = nodeMap.find(wedges[i]);
			wedgeProb1 = baseListProbCalc(u, w, nIt->second.neighbours);
			nIt = nodeMap.find(wedges[j]);
			wedgeProb2 = baseListProbCalc(u, w, nIt->second.neighbours);
			if (wedgeProb1 * wedgeProb2 >= t) {
				count++;
			}
		}
	}
	return count;
}

/***
* Greater Than Binary Search (first element which has value greater than val)
* wedges: Wedge Vector
* top, bottom: top and bottom index of the search
* val: target value
*/
int gtBinarySearch(vector<double> &wedges, int top, int bottom, double val) {
	int mid = (bottom - top)/2 + top;
	
	if (top > bottom) {
		return -1;
	}
	if (top == bottom) { //Squeeze
		if (wedges[mid] >= val) {
			if (wedges[mid + 1] < val) {
				return mid;
			}
		}
		return -1;
	}
	
	if (wedges[mid] >= val) {
		if (wedges[mid + 1] < val) {
			return mid;
		} else {
			return gtBinarySearch(wedges, mid + 1, bottom, val);
		}
	} else {
		return gtBinarySearch(wedges, top, mid - 1, val);
	}
}

/***
* Improved method of counting the number of uncertain butterflies from a list
*/
long imprListCount(vector<double> &wedges) {
	
	//It's cheaper to sort here than the keep a sorted list (presumably)
	sort(wedges.begin(), wedges.end(), [](const double a, const double b) {return a > b; });
	
	int i = 0;
	int j = 1;
	long count = 0;
	int n = wedges.size();
	bool flag = false;
	
	if (n < 2) {
		return 0;
	}
	
	for (; j < n; j++) {
		
		if (wedges[i] * wedges[j] < t) {
			flag = true;
			
			double minT = t / wedges[j]; //Minimum i wedge prob to satisfy t with j wedge
			
			i = gtBinarySearch(wedges, 0, i - 1, minT);
			
			if (i == -1) { //Terminate
				return count;
			}
			
		}
		
		count += (long) (i + 1);
		if (!flag) {
			i++;
		}
	}
	return count;
}

/***
* The UBFC Algoritm (our baseline solution)
*/
long ubfc() {
	long ubfCount = 0;
	unordered_map<int, Node>::iterator uIt, vIt;
	unordered_map<int, vector<int> > wedgeMap; //Stores wedge (via the middle node), the key is the end node
	unordered_map<int, vector<int> >::iterator wMapIt;
	
	vector<int> wTemp;
	for (int u : vID) {
		wedgeMap.clear();
		uIt = nodeMap.find(u);
		//cout << u << "\n";
		for (pair<int, double> vP : uIt->second.neighbours) {
			int v = vP.first;
			vIt = nodeMap.find(v);
			//cout << "\t" << v << "\n";
			if (vpComp(uIt->second, vIt->second)) {
				for (pair<int, double> wP : vIt->second.neighbours) {
					int w = wP.first;
					if (u == w) {
						continue;
					}
					//cout << "\t\t" << w << "\n";
					if (vpComp(uIt->second, nodeMap.find(w)->second)) {
						wMapIt = wedgeMap.find(w);
						if (wMapIt == wedgeMap.end()) {
							wTemp.clear();
							wTemp.push_back(v);
							wedgeMap.insert(make_pair(w, wTemp));
						} else {
							wMapIt->second.push_back(v);
						}
					} else {
						break;
					}
				}
			} else {
				break;
			}
		}
		for (pair<int, vector<int> > hashPair : wedgeMap) {
			ubfCount += baseListCount(u, hashPair.first, hashPair.second);
		}
	}
	return ubfCount;
}

/***
* Inserts a wedge (double) version
*/
void insertWedge(unordered_map<int, vector<double> > &wMap, int target, double prob) {
	unordered_map<int, vector<double> >::iterator wMapIt = wMap.find(target);
	if (wMapIt == wMap.end()) {
		vector<double> wTemp;
		wTemp.push_back(prob);
		wMap.insert(make_pair(target, wTemp));
	} else {
		wMapIt->second.push_back(prob);
	}
}

/***
* Counts the number of wedges in a single wedge map (set start-node)
*/
long wMapCount(unordered_map<int, vector<double> > &wMap) {
	long count = 0;
	for (pair<int, vector<double> > hashPair : wMap) {
		count += imprListCount(hashPair.second);
	}
	return count;
}

/***
* IUBFC-Vertex Priority first
*/
int iubfcvp() {
	long ubfCount = 0;
	unordered_map<int, Node>::iterator uIt, vIt;
	unordered_map<int, vector<double> > wedgeMap; //Stores wedge via ext prob
	unordered_map<int, vector<double> >::iterator wMapIt;
	double uvProb, vwProb, wProb; //Probs for e(u, v), e(v, w), w(u, v, w)
	double pruneT = t / (maxT * maxT); //Effective Pruning Threshold
	
	for (int u : vID) {
		wedgeMap.clear();
		uIt = nodeMap.find(u);
		for (pair<int, double> vP : uIt->second.neighbours) {
			int v = vP.first;
			vIt = nodeMap.find(v);
			if (vpComp(uIt->second, vIt->second)) {
				uvProb = vP.second;
				
				for (pair<int, double> wP : vIt->second.neighbours) {
					int w = wP.first;
					if (u == w) {
						continue;
					}
					if (vpComp(uIt->second, nodeMap.find(w)->second)) {
						vwProb = wP.second;
						if ((wProb = uvProb * vwProb) >= pruneT) {
							insertWedge(wedgeMap, w, wProb);
						}
						
					} else {
						break;
					}
				}
			} else {
				break;
			}
		}
		ubfCount += wMapCount(wedgeMap);
	}
	
	return ubfCount;
}

/***
* IUBFC-Existential Probability first
*/
long iubfcep() {
	long ubfCount = 0;
	unordered_map<int, Node>::iterator uIt, vIt;
	unordered_map<int, vector<double> > wedgeMap; //Stores wedge via ext prob
	unordered_map<int, vector<double> >::iterator wMapIt;
	double uvProb, vwProb, wProb; //Probs for e(u, v), e(v, w), w(u, v, w)
	double pruneT = t / (maxT * maxT); //Effective Pruning Threshold
	
	for (int u : vID) {
		wedgeMap.clear();
		uIt = nodeMap.find(u);
		for (pair<int, double> vP : uIt->second.neighbours) {
			int v = vP.first;
			vIt = nodeMap.find(v);
			uvProb = vP.second;
			if (uvProb >= pruneT) {
				for (pair<int, double> wP : vIt->second.neighbours) {
					int w = wP.first;
					if (u == w) {
						continue;
					}
					vwProb = wP.second;
					wProb = uvProb * vwProb;
					if (wProb >= pruneT) {
						if (vpComp(uIt->second, vIt->second) && vpComp(uIt->second, nodeMap.find(w)->second)) {
							insertWedge(wedgeMap, w, wProb);
						}
					} else {
						break;
					}
				}
			} else {
				break;
			}
		}
		ubfCount += wMapCount(wedgeMap);
	}
	
	return ubfCount;
}

/***
* Uncertain Butterfly Local Search - Node
*/
long ublsv(int node) {
	Node u = nodeMap.find(node)->second;
	unordered_map<int, vector<double> > wedgeMap;
	unordered_map<int, vector<double> >::iterator wMapIt;
	vector<double> wTemp;
	double uvProb, wProb;
	unordered_map<int, Node>::iterator uIt, vIt;
	
	for (pair<int, double> vP : u.neighbours) {
		int v = vP.first;
		vIt = nodeMap.find(v);
		uvProb = vP.second;
		if (uvProb < t) {
			continue;
		}
		for (pair<int, double> wP : vIt->second.neighbours) {
			int w = wP.first;
			if (node == w) {
				continue;
			}
			if ((wProb = uvProb * wP.second) >= t) {
				wMapIt = wedgeMap.find(w);
				if (wMapIt == wedgeMap.end()) {
					wTemp.clear();
					wTemp.push_back(wProb);
					wedgeMap.insert(make_pair(w, wTemp));
				} else {
					wMapIt->second.push_back(wProb);
				}
			}
		}
	}
	return wMapCount(wedgeMap);
}


/***
* Intersection of two sorted vectors (certain variation)
* Sort beforehand if neighbour lists are unsorted
*/
long certainVecIntersection(vector<pair<int, double> > &vec1, vector<pair<int, double> > &vec2, int pivot) {
	vector<pair<int, double> >::iterator v1It = vec1.begin();
	vector<pair<int, double> >::iterator v2It = vec2.begin();
	long count = 0;
	while (v1It != vec1.end() && v2It != vec2.end()) {
		if ((*v1It).first == pivot) {
			v1It++;
			continue;
		}
		if ((*v2It).first == pivot) {
			v2It++;
			continue;
		}
		if ((*v1It).first == (*v2It).first) {
			count ++;
			v1It++;
			v2It++;
		} else if ((*v1It).first < (*v2It).first) {
			v1It++;
		} else {
			v2It++;
		}
	}
	return count;
}

/***
* Intersection of two sorted vectors (uncertain variation)
* Sort beforehand if neighbour lists are unsorted
*/
long vecIntersection(vector<pair<int, double> > &vec1, vector<pair<int, double> > &vec2, int pivot, double exWProb) {
	vector<pair<int, double> >::iterator v1It = vec1.begin();
	vector<pair<int, double> >::iterator v2It = vec2.begin();
	long count = 0;
	while (v1It != vec1.end() && v2It != vec2.end()) {
		if ((*v1It).first == pivot) {
			v1It++;
			continue;
		}
		if ((*v2It).first == pivot) {
			v2It++;
			continue;
		}
		if ((*v1It).first == (*v2It).first) {
			if (((*v1It).second * (*v2It).second * exWProb) >= t) {
				count ++;
			}
			v1It++;
			v2It++;
		} else if ((*v1It).first < (*v2It).first) {
			v1It++;
		} else {
			v2It++;
		}
	}
	return count;
}

/***
* Uncertain Butterfly Local Search - Edge
*/
long ublse(int edge) {
	Edge e = allEdges[edge];
	if (e.prob < t) {
		return 0;
	}
	int n1 = e.n1;
	int n2 = e.n2;
	Node u = nodeMap.find(n1)->second;
	Node v = nodeMap.find(n2)->second;
	
	//Swap nodes so we enumerate over the smaller neighbour list
	if (u.neighbours.size() < v.neighbours.size()) {
		Node tempNode = v;
		v = u;
		u = tempNode;
	}
	
	long count = 0;
	for (pair<int, double> wP : v.neighbours) {
		if (wP.first != u.id) {
			Node w = nodeMap.find(wP.first)->second;
			double wProb = e.prob * wP.second;
			if (wProb < t) {
				continue;
			}
			count += vecIntersection(u.neighbours, w.neighbours, v.id, wProb);
		}
	}
	return count;
}

/***
* Factorial Funtion
*/
ll factorial(int n) {
	ll res = 1;
	for (ll i = 2; i <= n; i++) {
		res = res * i;
	}
	return res;
}

/***
* nCr function
*/
ll choose(int n, int r) {
	ll res = factorial(n) / (factorial(r) * factorial (n - r));
	return res;
}

ll lazyChoose(int n) {
	ll res = 0;
	for (int i = 1; i < n; i++) {
		res += (ll)(n - i);
	}
	return res;
}

ll scalableChoose(int n) {
	if (n == 1 || n == 0) {
		return 0;
	}
	if (n <= 17) {
		return choose(n, 2);
	} else {
		return lazyChoose(n);
	}
	
}

/***
* Calculates the number of certain butterflies with a start/end vertex
*/
ll wCertain(unordered_map<int, int> &map) {
	ll count = 0;
	unordered_map<int, int>::iterator mIt = map.begin();
	for (; mIt != map.end(); mIt++) {
		count += scalableChoose(mIt->second);
	}
	return count;
}

/***
* Certain vertex-centric butterfly enumeration for a given node
*/
ll blsv(int node) {
	unordered_map<int, Node>::iterator nodeIt = nodeMap.find(node);
	if (nodeIt == nodeMap.end()) {
		return 0;
	}
	Node u = nodeMap.find(node)->second;
	unordered_map<int, int> wedgeMap;
	unordered_map<int, int>::iterator wMapIt;
	unordered_map<int, Node>::iterator uIt, vIt;
	
	for (pair<int, double> vP : u.neighbours) {
		int v = vP.first;
		vIt = nodeMap.find(v);
		for (pair<int, double> wP : vIt->second.neighbours) {
			int w = wP.first;
			if (node == w) {
				continue;
			}
			wMapIt = wedgeMap.find(w);
			if (wMapIt == wedgeMap.end()) {
				wedgeMap.insert(make_pair(w, 1));
			} else {
				wMapIt->second++;
			}
		}
	}
	return wCertain(wedgeMap);
}

/***
* Certain edge-centric butterfly enumeration
*/
long blse(int edge) {
	Edge e = allEdges[edge];
	int n1 = e.n1;
	int n2 = e.n2;
	Node u = nodeMap.find(n1)->second;
	Node v = nodeMap.find(n2)->second;
	
	//Swap nodes so we enumerate over the smaller neighbour list
	if (u.neighbours.size() < v.neighbours.size()) {
		Node tempNode = v;
		v = u;
		u = tempNode;
	}
	
	long count = 0;
	for (pair<int, double> wP : v.neighbours) {
		if (wP.first != u.id) {
			Node w = nodeMap.find(wP.first)->second;
			count += certainVecIntersection(u.neighbours, w.neighbours, v.id);
		}
	}
	return count;
}


/***
* Returns a non-previously sampled random int within a range
*/
int rSample(int size, unordered_map<int, bool> &map) {
	int randVal = rand() % size;
	unordered_map<int, bool>::iterator mIt = map.find(randVal);
	if (mIt == map.end()) {
		map.insert(make_pair(randVal, true));
		return randVal;
	} else {
		return rSample(size, map);	
	}
}

/***
* Returns a random int within a range
*/
int rSample(int size) {
	return rand() % size;
}

/***
* Returns a random int wihtin a range excluding a specific number
*/
int rSampleXN(int size, int exc) {
	int randVal = rand() % size;
	if (randVal == exc) {
		return rSampleXN(randVal, exc);
	} else {
		return randVal;
	}
}

/***
* Uncertain Butterfly Sampling Algorithm - Vertex
*/
ld ubsv(int sampleNum) {
	unordered_map<int, bool> sampledMap;
	ld eUBFC = 0.0;
	for (int i = 0; i < sampleNum; i++) {
		int sampNode = rSample(vID.size(), sampledMap);
		ld localCount = (ld) ublsv(vID[sampNode]);
		ld extLocalCount = (localCount * (ld) vID.size())/4.0;
		eUBFC = (((ld)i * eUBFC) + extLocalCount) / (ld)(i + 1);
	}
	return eUBFC;
}

/***
* Uncertain Butterfly Sampling Algorithm - Edge
*/
ld ubse(int sampleNum) {
	unordered_map<int, bool> sampledMap;
	ld eUBFC = 0.0;
	for (int i = 0; i < sampleNum; i++) {
		int sampEdge = rSample(allEdges.size(), sampledMap);
		ld localCount = (ld) ublse(sampEdge);
		ld extLocalCount = localCount * (ld) allEdges.size() / 4.0;
		eUBFC = (((ld)i * eUBFC) + extLocalCount) / (ld)(i + 1);
	}
	return eUBFC;
}

/***
* Calcs the total number of certain wedges in the graph
*/
ll calcTotalWedges() {
	ll count = 0;
	for (auto &a : nodeMap) {
		count += scalableChoose(a.second.neighbours.size());
	}
	return count;
}

/***
* Proportion Estimation Method for Vertices
*/
ld pesv(int sampleNum, double alpha) {
	unordered_map<int, bool> sampledMap;
	ld eBFC = 0.0;
	for (int i = 0; i < sampleNum; i++) {
		int sampNode = rSample(vID.size(), sampledMap);
		ld localCount = (ld) blsv(sampNode);
		ld extLocalCount = (localCount * (ld) vID.size())/4.0;
		eBFC = (((ld)i * eBFC) + extLocalCount) / (ld)(i + 1);
	}
	return (ld)alpha * eBFC;
}

/***
* Proportion Estimation Method for Edges
*/
ld pese(int sampleNum, double alpha) {
	unordered_map<int, bool> sampledMap;
	ld eBFC = 0.0;
	for (int i = 0; i < sampleNum; i++) {
		int sampEdge = rSample(allEdges.size(), sampledMap);
		ld localCount = (ld) blse(sampEdge);
		ld extLocalCount = (localCount * (ld) allEdges.size())/4.0;
		eBFC = (((ld)i * eBFC) + extLocalCount) / (ld)(i + 1);
	}
	return (ld)alpha * eBFC;
}

/***
* Estimates the proportion based on the population sampling
*/
double propEs() {
	ld x = 2.56 * 2.56 * 0.25 / 0.0001 ;
	ld edgeC = (ld) allEdges.size();
	ld pop = edgeC * edgeC * edgeC * edgeC;
	ld n = pop * x / (pop + x - 1.0);
	
	int sampleSize = (int) n;
	//int sampleSize = 1000000;
	int suc = 0;
	unordered_map<int, bool> checkMap;
	for (int i = 0; i < sampleSize; i++) {
		checkMap.clear();
		double propE = allEdges[rSample(allEdges.size(), checkMap)].prob;
		if (propE < t) {
			continue;
		}
		propE = propE * allEdges[rSample(allEdges.size(), checkMap)].prob;
		if (propE < t) {
			continue;
		}
		propE = propE * allEdges[rSample(allEdges.size(), checkMap)].prob;
		if (propE < t) {
			continue;
		}
		propE = propE * allEdges[rSample(allEdges.size(), checkMap)].prob;
		if (propE < t) {
			continue;
		}
		suc++;
	}
//	cout << "Prop Sample Size: " << sampleSize << "\n";
//	cout << "Successes: " << suc << "\n";
	return (double) suc / (double) sampleSize;
}

/***
* argv[0] = ubfc
* argv[1] = Process Type
*	0 : UBFC (Baseline)
*	1 : IUBFC (Improved) - VP
*	2 : IUBFC (Improved) - EP
*	3 : UBS-Vertex
*	4 : UBS-Edge
*	5 : PES-Vertex
*	6 : PES-Edge
*	13 : Multiple UBS-Vertex Runs
*	14 : Multiple UBS-Edge Runs
*	15 : Multiple PES-Vertex Runs
*	16 : Multiple PES-Edge Runs
* argv[2] = Threshold Probability
* argv[3] = ID File
* argv[4] = Edge File
* argv[5] = # of Samples (optional)
* argv[6] = Output File (optional)
*/
int main(int argc, char *argv[]) {
	if (argc != 5 && argc != 6 && argc != 7) {
		printf("Insufficient Args\n");
		exit(1);
	}
	int proc = stoi(argv[1]);
	t = stod(argv[2]);
	if (t < 0.0 || t > 1.0) {
		printf("Incorrect Threshold\n");
		exit(1);
	}
	srand((unsigned int)time(NULL));
	auto start = chrono::high_resolution_clock().now();
	
	initialiseID(argv[3]);
	
	if (proc == 0) { //UBFC (BASELINE)
		cout << "BASELINE\n";
		initialiseEdge(argv[4]);
		initialiseVP();	
		cout << "UBF Count: " << ubfc() << "\n";
	} else if (proc == 1 || proc == 2) {
		initialiseEdgePrune(argv[4]);
		if (proc == 1) { //IUBFC-VP
			cout << "IUBFC-VP\n";
			initialiseVP();
			cout << "UBF Count: " << iubfcvp() << "\n";
		} else { //IUBFC-EP
			cout << "IUBFC-EP\n";
			initialiseEP();
			cout << "UBF Count: " << iubfcep() << "\n";
		}
	} else if (proc == 3 || proc == 4 || proc == 13 || proc == 14) { //UBS Approach
		int sampleSize;
		if (argc == 5) {
			sampleSize = 1000;
		} else {
			sampleSize = stoi(argv[5]);
		}
		initialiseEdge(argv[4]);
		if (proc == 3) {
			cout << "UBS-V\n";
			auto neststart = chrono::high_resolution_clock().now();
			cout << "UBF Estimation: " << ubsv(sampleSize) << "\n";
			auto neststop = chrono::high_resolution_clock().now();
			auto nestduration = chrono::duration_cast<chrono::milliseconds>(neststop - neststart).count();
			cout << "Sample Time: " << nestduration << "\n";
		} else if (proc == 4) {
			cout << "UBS-E\n";
			auto neststart = chrono::high_resolution_clock().now();
			cout << "UBF Estimation: " << ubse(sampleSize) << "\n";
			auto neststop = chrono::high_resolution_clock().now();
			auto nestduration = chrono::duration_cast<chrono::milliseconds>(neststop - neststart).count();
			cout << "Sample Time: " << nestduration << "\n";
		} else if (proc == 13) {
			for (int i = 1; i <= 500; i++) {
				cout << "Sample: " << i << "\n";
				auto neststart = chrono::high_resolution_clock().now();
				ld nestcount = ubsv(sampleSize);
				auto neststop = chrono::high_resolution_clock().now();
				auto nestduration = chrono::duration_cast<chrono::milliseconds>(neststop - neststart).count();
				string writeString = to_string(i);
				writeString += ",";
				writeString += to_string(nestcount);
				writeString += ",";
				writeString += to_string(nestduration);
				writeToFile(writeString, argv[6]);
			}
		} else if (proc == 14) {
			for (int i = 1; i <= 500; i++) {
				cout << "Sample: " << i << "\n";
				auto neststart = chrono::high_resolution_clock().now();
				ld nestcount = ubse(sampleSize);
				auto neststop = chrono::high_resolution_clock().now();
				auto nestduration = chrono::duration_cast<chrono::milliseconds>(neststop - neststart).count();
				string writeString = to_string(i);
				writeString += ",";
				writeString += to_string(nestcount);
				writeString += ",";
				writeString += to_string(nestduration);
				writeToFile(writeString, argv[6]);
			}
		}
	} else if (proc == 5 || proc == 6 || proc == 15 || proc == 16) { // PES Approach
		int sampleSize;
		if (argc == 5) {
			sampleSize = 1000;
		} else {
			sampleSize = stoi(argv[5]);
		}
		initialiseEdge(argv[4]);
		double alpha = propEs();
		cout << "Alpha: " << alpha << "\n";
		if (proc == 5) {
			cout << "PES-V\n";
			auto neststart = chrono::high_resolution_clock().now();
			cout << "UBF Estimation: " << pesv(sampleSize, alpha) << "\n";
			auto neststop = chrono::high_resolution_clock().now();
			auto nestduration = chrono::duration_cast<chrono::milliseconds>(neststop - neststart).count();
			cout << "Sample Time: " << nestduration << "\n";
		} else if (proc == 6) {
			cout << "PES-E\n";
			auto neststart = chrono::high_resolution_clock().now();
			cout << "UBF Estimation: " << pese(sampleSize, alpha) << "\n";
			auto neststop = chrono::high_resolution_clock().now();
			auto nestduration = chrono::duration_cast<chrono::milliseconds>(neststop - neststart).count();
			cout << "Sample Time: " << nestduration << "\n";
		} else if (proc == 15) {
			for (int i = 1; i <= 500; i++) {
				cout << "Sample: " << i << "\n";
				auto neststart = chrono::high_resolution_clock().now();
				double nestalpha = propEs();
				ld nestcount = pesv(sampleSize, nestalpha);
				auto neststop = chrono::high_resolution_clock().now();
				auto nestduration = chrono::duration_cast<chrono::milliseconds>(neststop - neststart).count();
				string writeString = to_string(i);
				writeString += ",";
				writeString += to_string(nestalpha);
				writeString += ",";
				writeString += to_string(nestcount);
				writeString += ",";
				writeString += to_string(nestduration);
				writeToFile(writeString, argv[6]);
			}
		} else if (proc == 16) {
			for (int i = 1; i <= 500; i++) {
				cout << "Sample: " << i << "\n";
				auto neststart = chrono::high_resolution_clock().now();
				double nestalpha = propEs();
				ld nestcount = pese(sampleSize, nestalpha);
				auto neststop = chrono::high_resolution_clock().now();
				auto nestduration = chrono::duration_cast<chrono::milliseconds>(neststop - neststart).count();
				string writeString = to_string(i);
				writeString += ",";
				writeString += to_string(nestalpha);
				writeString += ",";
				writeString += to_string(nestcount);
				writeString += ",";
				writeString += to_string(nestduration);
				writeToFile(writeString, argv[6]);
			}
		}
	} else {
		printf("Wrong Proc\n");
		exit(1);
	}
	auto stop = chrono::high_resolution_clock().now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start).count();
	cout << "Total Runtime: " << duration << "\n";
	
	if (argc == 7) {
		string timeString = "Total Runtime: ";
		timeString += to_string(duration);
		writeToFile(timeString, argv[6]);
	}
	
	return 0;
}
