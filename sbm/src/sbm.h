// functions used by both sbm and sbmt

using namespace Rcpp;

//-----------------------------------Global vars-------------------------------------

// Not modified
extern const int SelfTwice; ///Make 2 if you want to count self-edges twice (convention) in the directed case

// Network-Dependent
extern std::vector<double> Degree;  // Degree of nodes in the network or in Degree if directed
extern std::vector<double> outDegree;
extern std::vector<int> Count;  // Degree of nodes in the network or in Degree if directed
extern std::vector<int> outCount;
extern std::vector<int> LastEmpty;
extern std::vector<int> outLastEmpty;

extern std::vector<int> CurrentState;
extern std::vector<int> BestState;
extern std::vector<int> TrueState; // This records the true communities if they exist read in from file

extern std::vector<int> BestCommVertices; //[MaxComms]
extern std::vector<double>BestCommStubs; //if directed, keeps tally of edges originating in a class
extern std::vector<double>BestCommEnds; //for directed, keeps tally of edges ending in a class
extern std::vector<double>BestCommStubsTotal; // For networks with time slices, total stubs over all time periods
extern std::vector<double>BestCommEndsTotal; //    For use with DegreeCorrect 2

extern std::vector<int> CurrentCommVertices; //[MaxComms]
extern std::vector<double> CurrentCommStubs;
extern std::vector<double> CurrentCommEnds;
extern std::vector<double>CurrentCommStubsTotal;
extern std::vector<double>CurrentCommEndsTotal;

extern std::vector<double> NeighborSet;  // lists the weight of edges to that comm; is vertex specific
extern std::vector<double> outNeighborSet;

extern std::vector<double> BestEdgeMatrix;
extern std::vector<double> CurrentEdgeMatrix;

extern std::vector<double> SelfEdgeCounter;

extern int Nodes, MaxComms, KLPerNetwork, DegreeCorrect, T;
extern bool Directed;
extern double MaxScore; // records the *weight* of self-edges (not doubled) so they can be counted correctly.

//------------------------------------------------------------------------------------------------------------------------------------------------

void Initialize(IntegerMatrix AdjList, NumericMatrix AdjListWeight); // initializes the data structures for KL

double ComputeInitialScore(); // computes the initial score after initialization

void RunKL(IntegerMatrix AdjList, NumericMatrix AdjListWeight, IntegerMatrix outAdjList, NumericMatrix outAdjListWeight); // runs Kernighan-Lin once.

void ComputeNeighborSet(int vertex, int option, IntegerMatrix AdjList, NumericMatrix AdjListWeight, IntegerMatrix outAdjList, NumericMatrix outAdjListWedght, int time);  // computes the neighbor set for a vertex, using either best or currentstates

double ComputeProposal(int vertex, int from, int destination); // computes the value of the particular change

void UpdateMatrices(int vertex, int option, int from, int destination); // this updates either the best
IntegerVector randomComms(int Nodes, int MaxComms); // generates random commmunity assignments

double LogFunction(double x); // this returns x*log(x) and zero if x=0
