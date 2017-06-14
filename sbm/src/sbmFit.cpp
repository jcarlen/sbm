/// TO DO: 1) Add twice diag option for undirected? Move to R side?
///        2) Fix the code for reading in given initial communities (on the R and cpp ends)
///        3) Allow network input and function converts to edgelist

/// CHANGES: 1) Edge matrices are vectors (get same behavior using row and column indexing) so they can be declared and then size changed.
///          2) Use Rcpp for random community generation, no rng dependence
///          3) Doubles instead of floats - need precision for LogFunction

// -------------------------------------------------------------

/// Jane Carlen
///
/// March, 18 2017
///
/// DIRECTED and undirected, degree-corrected stochastic block model implementation
///
/// Adapted from Brian Karrer http://www-personal.umich.edu/~mejn/dcsbm/KLOptimization.cpp
///
/// Real-valued, non-neg weights accepted via weights arg (rather than repeated edges in list)
///
/// This version is similar to KLOptimization_directed, but connects to R via RCPP


// -------------------------------------------------------------

// Please do not distribute without contacting karrerb@umich.edu.

// This code uses the GNU GSL random number generators.  If these are unavailable, please
// replace the random number generator.  The default is to seed the random number generator using
// the system clock and to use the mersenne twister.

// Please read the implementation directions below to start using the program.

/*#include <stdio.h>
#include <string.h>

#include <fstream>
#include <limits>


#include <ctime>*/


//#include <math.h>
//#include <vector>
#include <Rcpp.h>
//#include "sbmFit.h" - was a file for randomComms function so I didn't have copy of it in thi sfile and KL.cpp

using namespace Rcpp;


//*********************** GLOBAL VARIABLES *******************************************************************

// Not modified
const int SelfTwice = 1; ///Make 2 if you want to count self-edges twice (convention) in the directed case

// Network-Dependent
std::vector<double> Degree;  // Degree of nodes in the network or in Degree if directed
std::vector<double> outDegree;
std::vector<int> Count;  // Degree of nodes in the network or in Degree if directed
std::vector<int> outCount;
std::vector<int> LastEmpty;
std::vector<int> outLastEmpty;

std::vector<int> CurrentState;
std::vector<int> BestState;
std::vector<int> TrueState; // This records the true communities if they exist read in from file

std::vector<int> BestCommVertices; //[MaxComms]
std::vector<double>BestCommStubs; //if directed, keeps tally of edges originating in a class
std::vector<double>BestCommEnds; //for directed, keeps tally of edges ending in a class

std::vector<int> CurrentCommVertices; //[MaxComms]
std::vector<double> CurrentCommStubs;
std::vector<double> CurrentCommEnds;

std::vector<int> NeighborIndex; // [MaxComms] //lists the comm
std::vector<int> outNeighborIndex;
std::vector<int> TempNeighborIndex;
std::vector<int> outTempNeighborIndex;

std::vector<double> NeighborSet;  // lists the weight of edges to that comm
std::vector<double> outNeighborSet;
std::vector<double> TempNeighborSet;
std::vector<double> outTempNeighborSet;

std::vector<double> BestEdgeMatrix;
std::vector<double> CurrentEdgeMatrix;

std::vector<int> ActualDiffComms;
std::vector<int> outActualDiffComms;
std::vector<double> SelfEdgeCounter;

int Nodes, MaxComms, KLPerNetwork, DegreeCorrect, T;
bool Directed;
double MaxScore; // records the *weight* of self-edges (not doubled) so they can be counted correctly.

//*********************** FUNCTION DECLARATIONS **************************************************************

void Initialize(IntegerMatrix AdjList, NumericMatrix AdjListWeight); // initializes the data structures for KL
void RunKL(IntegerMatrix AdjList, NumericMatrix AdjListWeight, IntegerMatrix outAdjList, NumericMatrix outAdjListWeight); // runs Kernighan-Lin once.
double ComputeInitialScore(); // computes the initial score after initialization
void ComputeNeighborSet(int vertex, int option, IntegerMatrix AdjList, NumericMatrix AdjListWeight, IntegerMatrix outAdjList, NumericMatrix outAdjListWeight, int time);  // computes the neighbor set for a vertex, using either best or currentstates
double ComputeProposal(int vertex, int from, int destination); // computes the value of the particular change
void UpdateMatrices(int vertex, int option, int from, int destination); // this updates either the best
double LogFunction(double x); // this returns x*log(x) and zero if x=0
IntegerVector randomComms(int Nodes, int MaxComms); // generates random commmunity assignments

// for calling KL directly, for time-dependent models
void Setup(int Nodes, int MaxComms, bool Directed);
void InitializeKLt(List AdjList, List AdjListWeight); // initializes the data structures for KL
List RunKLt(SEXP edgelistTime, const IntegerMatrix & edgelistTotal, const NumericVector & weightsTotal, const int maxComms, const bool directed, const int klPerNetwork, const bool degreeCorrect, const int nodes);

//*********************** MAIN PROGRAM ***********************************************************************

// Main program

// [[Rcpp::export]]

List sbmFit(const IntegerMatrix & edgelist, const int maxComms, const int degreeCorrect, const bool directed, const int klPerNetwork, const NumericVector weights) {  //, const IntegerVector seedComms)

    int i, j, k;
    
    // Get number of nodes (and create NodeSet, but might not need it)
    std::set<int> NodeSet;
    int counter = edgelist.nrow();
    
    for (i = 1; i < counter; i++)
    {
      NodeSet.insert(edgelist(i, 0));
      NodeSet.insert(edgelist(i, 1));
    }
    
    Nodes = NodeSet.size();
    std::cout << "Nodes " << Nodes << std::endl;
    MaxComms = maxComms;
    DegreeCorrect = degreeCorrect;
    Directed = directed;
    KLPerNetwork = klPerNetwork;
    T = 1; //time steps, if applicable
    
    /// Make space in Global Vars
    CurrentState.assign(Nodes, 0);
    BestState.assign(Nodes, 0);
    TrueState.assign(Nodes, 0); // This records the true communities if they exist read in from the file
    
    Degree.assign(Nodes, 0.0);  // Degree of nodes in the network or in Degree if directed
    Count.assign(Nodes, 0);  // Degree of nodes in the network or in Degree if directed
    LastEmpty.assign(Nodes, 0);
    
    BestCommVertices.assign(MaxComms, 0); //[MaxComms]
    BestCommStubs.assign(MaxComms, 0.0); //if directed, keeps tally of edges originating in a class
    BestEdgeMatrix.assign(MaxComms*MaxComms, 0.0);
    
    CurrentCommVertices.assign(MaxComms, 0); //[MaxComms]
    CurrentCommStubs.assign(MaxComms, 0.0); //if directed, keeps tally of edges originating in a class
    CurrentEdgeMatrix.assign(MaxComms*MaxComms, 0.0);
    
    NeighborIndex.assign(MaxComms, 0); // [MaxComms] //lists the comm
    TempNeighborIndex.assign(MaxComms, 0);
    NeighborSet.assign(MaxComms, 0.0);  // lists the weight of edges to that comm
    TempNeighborSet.assign(MaxComms, 0.0);
    
    ActualDiffComms.assign(T, 0);
    
    // Set up network structure vars - Degree, Count, AdjLists
    if (!Directed)
    {
        
        for(i=0; i < counter; i++) // First we count the degrees by scanning through the list once
        {
            Degree[edgelist(i, 0)]+=weights[i];
            Degree[edgelist(i, 1)]+=weights[i];
            Count[edgelist(i, 0)]++;
            Count[edgelist(i, 1)]++;
            
        }
    }
 
    int outmaxCount = 1; //have to set a default for directed.
    
    if (Directed)
    {
        
        // Directed-only Initializations
        outDegree.assign(Nodes, 0.0);
        outCount.assign(Nodes, 0);
        outLastEmpty.assign(Nodes, 0);
        
        BestCommEnds.assign(MaxComms, 0.0); //for directed, keeps tally of edges ending in a class
        CurrentCommEnds.assign(MaxComms, 0.0); //for directed, keeps tally of edges ending in a class
        
        outNeighborIndex.assign(MaxComms, 0);
        outTempNeighborIndex.assign(MaxComms, 0);
        outNeighborSet.assign(MaxComms, 0.0);
        outTempNeighborSet.assign(MaxComms, 0.0);
        
        outActualDiffComms.assign(T, 0);
        
        for(i=0; i < counter; i++)
        {
            Degree[edgelist(i, 1)]+=weights[i];
            outDegree[edgelist(i, 0)]+=weights[i];
            Count[edgelist(i, 1)]++;
            outCount[edgelist(i, 0)]++;
        }
        
        outmaxCount = *std::max_element(outCount.begin(), outCount.end());
    }

    // Create big enough adjlist matrices (less efficient than dynamic mem cpp-only implementation, but seems necessary)
    // keep at this scope so it's set correctly
    // undirected & directed
    int maxCount = *std::max_element(Count.begin(), Count.end());
    IntegerMatrix AdjList(Nodes, maxCount);
    NumericMatrix AdjListWeight(Nodes, maxCount);
    // directed only
    IntegerMatrix outAdjList(Nodes, outmaxCount);
    NumericMatrix outAdjListWeight(Nodes, outmaxCount);
    
    //Fill in adj matrices
    if (!Directed)
    {

        for(i=0; i < counter; i++)
        {
            
            AdjList(edgelist(i, 0), LastEmpty[edgelist(i, 0)]) = edgelist(i, 1);
            AdjListWeight(edgelist(i, 0), LastEmpty[edgelist(i, 0)]) = weights[i];
            LastEmpty[edgelist(i, 0)]+=1;
            
            AdjList(edgelist(i, 1), LastEmpty[edgelist(i, 1)]) = edgelist(i, 0);
            AdjListWeight(edgelist(i, 1), LastEmpty[edgelist(i, 1)]) = weights[i];
            LastEmpty[edgelist(i, 1)]+=1;
        }
        
        
    }
    
    if (Directed) {
        
        for(i=0; i < counter; i++)
        {
            
            AdjList(edgelist(i, 1), LastEmpty[edgelist(i, 1)]) = edgelist(i, 0);
            AdjListWeight(edgelist(i, 1), LastEmpty[edgelist(i, 1)]) = weights[i];
            LastEmpty[edgelist(i, 1)]+=1;
            
            outAdjList(edgelist(i, 0), outLastEmpty[edgelist(i, 0)]) = edgelist(i, 1);
            outAdjListWeight(edgelist(i, 0), outLastEmpty[edgelist(i, 0)]) = weights[i];
            
            outLastEmpty[edgelist(i, 0)]+=1;
        }
        
    }

    double HighestScore = -std::numeric_limits<double>::max( );
    /*double VIValue = 0;
    double NMIValue = 0;*/
    
    // For reporting best state
    std::vector<int> SavedState(Nodes, 0);
    std::vector<int> SavedCommVertices(MaxComms, 0);
    std::vector<double> SavedCommStubs(MaxComms, 0.0);
    std::vector<double> SavedCommEnds(MaxComms, 0.0);
    std::vector<double> SavedEdgeMatrix(MaxComms*MaxComms, 0.0);
    
    for(j=0; j < KLPerNetwork; j++)
    {
        
        RunKL(AdjList, AdjListWeight, outAdjList, outAdjListWeight);
        
        if(MaxScore >= HighestScore)
        {
            HighestScore = MaxScore;
            /*if(TrueCommsAvailable == 1)
             {
                VIValue = ComputeVI();
                NMIValue = ComputeNMI();
             }*/
            
            for(i=0; i < MaxComms; i++)
                
            {
                SavedCommVertices[i] = BestCommVertices[i];
                //Rcout << BestCommVertices[i] << ", ";
                SavedCommStubs[i] = BestCommStubs[i];
                if (Directed) SavedCommEnds[i] = BestCommEnds[i];
                    for(k=0; k < MaxComms; k++)
                        SavedEdgeMatrix[i*MaxComms+k] = BestEdgeMatrix[i*MaxComms+k];
            }
             for(i=0; i < Nodes; i++)
                 SavedState[i] = BestState[i];
         }
    }
    
    // because PrintResults are written for best values we copy them
    // back over from the saved values which are the best ones.
    for(i=0; i < MaxComms; i++)
    {
      BestCommVertices[i] = SavedCommVertices[i];
      BestCommStubs[i] = SavedCommStubs[i];
      if (Directed) SavedCommEnds[i] = BestCommEnds[i];
      for(k=0; k < MaxComms; k++)
          BestEdgeMatrix[i*MaxComms+k] = SavedEdgeMatrix[i*MaxComms+k];
    }
    
    for(i=0; i < Nodes; i++)
        BestState[i] = SavedState[i];
    
    return List::create(Rcpp::Named("FoundComms") = BestState, //note they will correspond to ORDER of IDs
                        Rcpp::Named("EdgeMatrix") = BestEdgeMatrix,
                        Rcpp::Named("HighestScore") = HighestScore);

 }


void RunKL(IntegerMatrix AdjList, NumericMatrix AdjListWeight, IntegerMatrix outAdjList, NumericMatrix outAdjListWeight)
{
    int i,j,k;
    int MaxIndex;
    double CurrentScore;  // records the current log-likelihood
    int MaxVertex;  // this records the index of the largest vertex ratio found so far
    double MaxRatio;  // records the value of the ratio, actually it's the log of the ratio
    int MaxPriority; // records the community that the vertex wants to go toI
    long int tempvertex = 0;
    std::vector<int> ChangeSet(Nodes, 0);
    std::vector<int> UpdateIndex(Nodes, -1);
    
    double prevMaxScore = -std::numeric_limits<double>::max( );
    long double tolerance = 0.00000001; // this prevents loops due to numerical errors.
    
    double ProposalRatio;
    double value;
    int Priority;
    
    Initialize(AdjList, AdjListWeight);
    
    // This returns the log of the initial score
    MaxScore = ComputeInitialScore();
    
    while(MaxScore >= prevMaxScore + tolerance)
    {
        Rcout << "MAX SCORE IS: " << MaxScore << std::endl;
        // we start with everything equal to the best values
        CurrentScore = MaxScore;
        prevMaxScore = MaxScore;
        MaxIndex = -1;
        
        // ChangeSet records which vertices are able to move, in that they haven't
        // already moved during this KL step.  Update index will tell when the vertex
        // was chosen to move.
        for(i=0; i < Nodes; i++)
        {
            CurrentState[i] = BestState[i];
            ChangeSet[i] = i;
            UpdateIndex[i] = -1;
        }

        for(i=0; i < MaxComms; i++)
        {
            CurrentCommVertices[i] = BestCommVertices[i];
            CurrentCommStubs[i] = BestCommStubs[i];
            if ( Directed ) CurrentCommEnds[i] = BestCommEnds[i];
            for(j=0; j < MaxComms; j++)
                CurrentEdgeMatrix[i*MaxComms+j] = BestEdgeMatrix[i*MaxComms+j];
        }
        
        // This loop moves each vertex once
        // Note that we DONT reinitialize changeset as this is unnecessary
        // This would make it a factor of 2 slower.
        for(i=0; i < Nodes; i++)
        {
            MaxVertex = 0;
            MaxRatio = -std::numeric_limits<double>::max( );
            MaxPriority = 0;
            
            // This loop selects which vertex to move
            for(j=0; j < Nodes-i; j++)
            {
                
                // get proposal and proposal ratio for ChangeSet[j]
                Priority = 0;
                ProposalRatio = -std::numeric_limits<double>::max( );
                // we first compute the neighbor set of the vertex, this is fixed
                // and the same for every change,
                // computing this first makes this more efficient
                // zero indicates run with current communities
                
                ComputeNeighborSet(ChangeSet[j], 0, AdjList, AdjListWeight, outAdjList, outAdjListWeight, 0); //time is 0. time steps is 1
                //Rcout << "Finished ComputeNeighborSet" << std::endl;
                
                for(k=0; k < MaxComms; k++)
                {
                    // we compute the value of vertex ChangeSet[j] going to k
                    // we DONT allow a vertex to remain where it was
                    // This is essential to enforce so that it will go downhill and not be greedy
                    if(k != CurrentState[ChangeSet[j]])
                    {
                        value = ComputeProposal(ChangeSet[j], CurrentState[ChangeSet[j]], k);
                        if(value > ProposalRatio)
                        {
                            Priority = k;
                            ProposalRatio = value;

                        }
                    }
                }
                
                //Rcout << "Finished ComputeProposal" << std::endl;
                
                // check whether its higher than what you already have as the max KL move
                if(ProposalRatio > MaxRatio)
                {
                    MaxVertex = j;  // Note this is not the vertex j, but the vertex given by ChangeSet[j]
                    MaxRatio = ProposalRatio;
                    MaxPriority = Priority;
                }
            }
            
            // now we move it, first recording the current neighbors so that
            // we can update the matrices properly
            ComputeNeighborSet(ChangeSet[MaxVertex], 0, AdjList, AdjListWeight, outAdjList, outAdjListWeight, 0);
            
            // This updates the matrices to represent the vertices new state
            UpdateMatrices(ChangeSet[MaxVertex], 0, CurrentState[ChangeSet[MaxVertex]], MaxPriority);
            
            CurrentState[ChangeSet[MaxVertex]] = MaxPriority;
            
            // we are using logs so we add the maxratio to the current score for the new score
            CurrentScore = CurrentScore + MaxRatio;
            
            UpdateIndex[ChangeSet[MaxVertex]] = i;
            // we switch it with the last element of changeset, removing it from further consideration
            // until we have moved the other vertices
            tempvertex = ChangeSet[MaxVertex];
            ChangeSet[MaxVertex] = ChangeSet[Nodes-i-1];
            ChangeSet[Nodes-i-1] = tempvertex;
            
            // now if the new state is better than the previous best state we record this
            // MaxIndex in combination with UpdateIndex
            // telling us where we are in the movement of vertices

            if(CurrentScore > MaxScore)
            {
                MaxScore = CurrentScore;
                MaxIndex = i;
            }
        }
        
        // now we update BestState if a change resulted in a higher maximum
        // by implementing all the changes found above
        
        // There is a potential for speeding this up here.
        if(MaxIndex != -1)
        {
            for(i=0; i < Nodes; i++)
            {
                // we only make the changes to beststate that happened before or equal to maxindex
                // no other vertex is updated
                // fortunately the update order is irrelevant to the final result so
                // we can just do it this way
                // Because we force all moves to be different, these updates are all a switch of community
                if(UpdateIndex[i] <= MaxIndex)
                {
                    // the option 1 does update corresponding to the best states / matrices
                    ComputeNeighborSet(i, 1, AdjList, AdjListWeight, outAdjList, outAdjListWeight, 0);
                    UpdateMatrices(i, 1, BestState[i], CurrentState[i]); // 1 does best matrix update
                    
                    BestState[i] = CurrentState[i];
                }
            }
        }
    }
    
    return;
}


void Initialize(IntegerMatrix AdjList, NumericMatrix AdjListWeight) // Calculate Vertices, Stubs, Ends, Matrix based on Comm initiation (currently random initiation only)
{
    int i, j;
    int neighbor;
    // double sum;
    IntegerVector randComms = randomComms(Nodes, MaxComms);
    
    if ( !Directed )
    {
        
        for(i=0; i < MaxComms; i++) ///initialize
        {
            BestCommVertices[i] = 0;
            BestCommStubs[i] = 0;
            for(j=0; j < MaxComms; j++)
                BestEdgeMatrix[i*MaxComms+j] = 0;
        }
       
        for(i=0; i < Nodes; i++)
        {
            BestState[i] = randComms[i];
            /*if(InitializationOption == 1)
                BestState[i] = TrueState[i];*/
            BestCommVertices[BestState[i]]++;
            BestCommStubs[BestState[i]] += Degree[i]; ///stubs will add to twice the total edge weight
        }
 
        /// sum = 0;
        
        for(i=0; i < Nodes; i++)
        {
            for(j=0; j < LastEmpty[i]; j++)
            {
                neighbor = AdjList(i, j); ///each edge listed twice, once for each end
                double neighborWeight = AdjListWeight(i, j);
                /// sum += neighborWeight;
                
                /// twice the diagonal of BestEdgeMatrix + everything else in the matrix once should equal = 2x edges
                
                if (BestState[neighbor] == BestState[i])
                    neighborWeight = (neighborWeight/2.0);
                
                BestEdgeMatrix[(BestState[i])*MaxComms + BestState[neighbor]] += neighborWeight;
                
            }
        }
    }
   
    
    if ( Directed )
    {
        for(i=0; i < MaxComms; i++) /// initialize
        {
            BestCommVertices[i] = 0;
            BestCommStubs[i] = 0;
            BestCommEnds[i] = 0;
            for(j=0; j < MaxComms; j++)
                BestEdgeMatrix[i*MaxComms+j] = 0;
        }
        
        for(i=0; i < Nodes; i++)
        {
        
            BestState[i] = randComms(i);
            /*if(InitializationOption == 1)
                BestState[i] = TrueState[i];*/
            BestCommVertices[BestState[i]]++; ///track number of vertices in each class
            BestCommStubs[BestState[i]] += outDegree[i];
            BestCommEnds[BestState[i]] += Degree[i];
        }
        
        // sum = 0;
        for(i=0; i < Nodes; i++)
        {
            for(j=0; j < LastEmpty[i]; j++)
            {
                neighbor = AdjList(i, j);
                
                BestEdgeMatrix[BestState[neighbor] * MaxComms + BestState[i]] += AdjListWeight(i, j);
                // sum += AdjListWeight(i, j); /// Note self-edges get counted ONCE because we're only interating over out-edges
                
                /// count self edges twice?
                if (i == neighbor && SelfTwice == 2) {
                    int position = BestState[i]*MaxComms + BestState[neighbor];
                    BestEdgeMatrix[position] += AdjListWeight(i, j);
                    // sum += AdjListWeight(i, j);
                }
                
                
                
            }
            
        }
    }
    
    /// also, need to double-check that sum is correct so everything adds to one
    
     return;
    
}


double ComputeInitialScore()
{
    // For the running of the KL algorithm itself this does not matter as all we use
    // are whether the score increases
    // We will want this when we compare different initializations
    
    int i,j;
    double sum = 0;
    
    if ( !Directed ) // this actually returns 1/2 the unnormalized log-likelihood listed in the paper
    {
        
        for (int t = 0; t < T; t++) {

            for(i=0; i < MaxComms; i++)
            {
                if(DegreeCorrect == 1)
                    sum = sum - LogFunction(BestCommStubs[i + t*MaxComms]); //LogFunction is x*log(x) or 0 if x is 0
                
                if(DegreeCorrect == 0)
                {
                    if(BestCommVertices[i] != 0)
                        sum = sum - (BestCommStubs[i + t*MaxComms]) * log(BestCommVertices[i]);
                }
                
                for(j=i; j < MaxComms; j++)
                {
                    if(j != i)
                        sum = sum + LogFunction(BestEdgeMatrix[i*MaxComms+j + t*MaxComms*MaxComms]);
                    
                    if(i==j)
                        sum = sum + .5*LogFunction(2*BestEdgeMatrix[i*MaxComms+j + t*MaxComms*MaxComms]);
                }
                
            }
        //Rcout << "sum after time" << t << " " << sum <<std::endl;
        }
        
    }
    
    if ( Directed )
    {
     
        for (int t = 0; t < T; t++) {
            
            for(i=0; i < MaxComms; i++)
            {
                for(j = 0; j < MaxComms; j++)
                {
                    double bestval;
                    
                    if (DegreeCorrect == 1) {
                        bestval = BestCommStubs[i + t*MaxComms] * BestCommEnds[j + t*MaxComms];
                    } else {
                        bestval = BestCommVertices[i] * BestCommVertices[j];
                    }
                    
                    sum = sum + LogFunction(BestEdgeMatrix[i*MaxComms+j + t*MaxComms*MaxComms]);
                    
                    if ( bestval > 0 && BestEdgeMatrix[i*MaxComms+j + t*MaxComms*MaxComms] > 0)
                        sum = sum - BestEdgeMatrix[i*MaxComms+j + t*MaxComms*MaxComms] * log(bestval);
                    
                }
                
            }
            
        }
    }
    
    Rcout << "sum " << sum << std::endl;
    return sum;
    
    
}

// We compute this using the current comm matrices
// We avoid the potential pitfalls of huge intermediate numbers by adding logs together.  So we treat 0 log 0 as 0.
// We return 0 for degree zero vertices (which really shouldn't be sent into the program
// in the first place.)
// We also return 0 for from = destination cause there is no change then.
// Here we use base e.  It returns the log of the actual value.
// Again this is half of the change in the unnormalized log-likelihood listed in the paper

void ComputeNeighborSet(int vertex, int option, IntegerMatrix AdjList, NumericMatrix AdjListWeight, IntegerMatrix outAdjList, NumericMatrix outAdjListWeight, int time)
{
    
    int i;
    int neighbor;
    double tol = std::numeric_limits<double>::epsilon();
    SelfEdgeCounter.assign(T, 0.0);
    
    if ( !Directed )
    {
        //Rcout << "ComputeNeighborSet_Undirected" << std::endl;
        
        for(i=0; i < MaxComms; i++)
        {
            TempNeighborIndex[i] = i;
            TempNeighborSet[i] = 0;
            NeighborIndex[i + time*MaxComms] = i;
            NeighborSet[i + time*MaxComms] = 0;
        }

        // NOTE SINCE A SELF-EDGE SHOWS UP TWICE IN THE ADJLIST THIS DOUBLE
        // COUNTS THESE EDGES, WE RECORD THE NUMBER OF TIMES THIS HAPPENS
        // IN A SEPARATE VARIABLE AND THEN DIVIDE BY TWO
        
        for(i=0; i < Count[vertex + time*Nodes]; i++)
        {
            
            neighbor = AdjList(vertex, i);
            
            if(neighbor != vertex)
            {

                if(option == 0)
                    TempNeighborSet[CurrentState[neighbor]] += AdjListWeight(vertex, i);
                
                if(option == 1)
                    TempNeighborSet[BestState[neighbor]] += AdjListWeight(vertex, i);
            }
            if(neighbor == vertex)
                SelfEdgeCounter[time] += AdjListWeight(vertex, i);
        }

        SelfEdgeCounter[time] = SelfEdgeCounter[time]/2; /// ONLY if undirected
        
        ActualDiffComms[time] = 0;
        for(i=0; i < MaxComms; i++)
        {
            if(TempNeighborSet[i] > tol)
            {
                NeighborIndex[ActualDiffComms[time] + time*MaxComms] = TempNeighborIndex[i];
                NeighborSet[ActualDiffComms[time] + time*MaxComms] = TempNeighborSet[i];
                ActualDiffComms[time]++;
            }
        }
        
        return;
    }
    
    if ( Directed )
    {
         //Rcout << "ComputeNeighborSet_Directed" << std::endl;
         for(i=0; i < MaxComms; i++)
         {
             TempNeighborIndex[i] = i;
             TempNeighborSet[i] = 0;
             outTempNeighborIndex[i] = i;
             outTempNeighborSet[i] = 0;
             NeighborIndex[i + time*MaxComms] = i;
             NeighborSet[i + time*MaxComms] = 0;
             outNeighborIndex[i + time*MaxComms] = i;
             outNeighborSet[i + time*MaxComms] = 0;
         }
     
     
         for(i=0; i < Count[vertex + time*Nodes]; i++)
         {
             neighbor = AdjList(vertex, i);
             if(neighbor != vertex)
             {
                 if(option == 0)
                     TempNeighborSet[CurrentState[neighbor]] += AdjListWeight(vertex, i);  /// the first entry lists the comm and the second entry lists the number of edges to that comm
         
                 if(option == 1)
                     TempNeighborSet[BestState[neighbor]] += AdjListWeight(vertex, i);
             }
             if(neighbor == vertex)
                 SelfEdgeCounter[time] += SelfTwice * AdjListWeight(vertex, i);  /// count self-edges ONCE unless selftwice is 2
         }
         
         for(i=0; i < outCount[vertex + time*Nodes]; i++)
         {
             neighbor = outAdjList(vertex, i);
             if(neighbor != vertex)
             {
                if(option == 0)
                    outTempNeighborSet[CurrentState[neighbor]] += outAdjListWeight(vertex, i);
     
                if(option == 1)
                    outTempNeighborSet[BestState[neighbor]] += outAdjListWeight(vertex, i);
             }
     
         }
     
         ActualDiffComms[time] = 0;
         for(i=0; i < MaxComms; i++)
         {
            if(TempNeighborSet[i] > tol )
            {
                NeighborIndex[ActualDiffComms[time] + time*MaxComms] = TempNeighborIndex[i];
                NeighborSet[ActualDiffComms[time] + time*MaxComms] = TempNeighborSet[i];
                ActualDiffComms[time]++;
            }
         }
        
         outActualDiffComms[time] = 0;
         for(i=0; i < MaxComms; i++)
         {
             if(outTempNeighborSet[i] > tol )
             {
                 outNeighborIndex[outActualDiffComms[time] + time*MaxComms] = outTempNeighborIndex[i];
                 outNeighborSet[outActualDiffComms[time] + time*MaxComms] = outTempNeighborSet[i];
                 outActualDiffComms[time]++;
             }
         }
     
     return;
     }
}

double ComputeProposal(int vertex, int from, int destination)
{
    int i;
    double ratiovalue = 0;
    double fromcount, destcount, outfromcount, outdestcount;
    
    double help1;
    double help2;
    double help3;
    double help4;
    
    if(from == destination)
        return 0;
    
    for (int t = 0; t < T; t ++)
    {
        fromcount = 0;
        destcount = 0;
        outfromcount = 0;
        outdestcount = 0;
        
         if (!Directed)
         {
             
             // if the degree of the vertex is zero we know nothing about it
             // in this case we don't ever change its community
             // at the end we put all degree zeroes into their own group
             if(DegreeCorrect == 1)
             {
                 if(Degree[vertex + t*Nodes] == 0)
                 return 0;
             }
             
             // we first add up all the cross-terms (between communities that are not from / destination)
             for(i=0; i < ActualDiffComms[t]; i++)
             {
                 // we lost NeighborSet[i] edges to NeighborIndex[i] from the from comm
                 // we gain the same amount in the destination comm
                 // IFF the comms were not from and destination
                 if( (NeighborIndex[i + t*MaxComms] != from) && (NeighborIndex[i + t*MaxComms] != destination))
                 {
                        // do update NOTE: each community mcc' gets updated once if it had edges switch out
                        // which is correct, remembering that mcc' is symmetric (///undirected case) and we only count c < c' here
             
                        help1 = double(CurrentEdgeMatrix[from * MaxComms + NeighborIndex[i + t*MaxComms] + t*MaxComms*MaxComms]);
                        help2 = double(CurrentEdgeMatrix[destination * MaxComms + NeighborIndex[i + t*MaxComms] + t*MaxComms*MaxComms]);
                        help3 = double(NeighborSet[i + t*MaxComms] + t*MaxComms);
                        if (help1 - help3 < 0) std::cout << "help1-help3  " << help1-help3 << std::endl;
                        if (help2 < 0) std::cout << "help2 " << help2 << std::endl;
                        if (help3 < 0) std::cout << "help3 " << help3 << std::endl;
                     
                     
                        ratiovalue = ratiovalue + LogFunction(help1-help3) - LogFunction(help1);
                        ratiovalue = ratiovalue + LogFunction(help2+help3) - LogFunction(help2);
                 }
             
                 if(NeighborIndex[i + t*MaxComms] == from)
                     fromcount = NeighborSet[i + t*MaxComms];
             
                 if(NeighborIndex[i + t*MaxComms] == destination)
                     destcount = NeighborSet[i + t*MaxComms];
             }
             
             // now we add in the term corresponding to from / dest
             
             help1 = double(CurrentEdgeMatrix[from * MaxComms + destination + t*MaxComms*MaxComms]);
             help2 = double(fromcount-destcount);
             if (help1 + help2 < 0) std::cout << "help1 + help2  " << help1 + help2 << std::endl;
             if (help1 < 0) std::cout << "help1 " << help1 << std::endl;
             
             ratiovalue = ratiovalue + LogFunction(help1 + help2) - LogFunction(help1);
             
             // now we add in the terms corresponding to from
             if(DegreeCorrect == 1)
             {
                 help1 = double(CurrentCommStubs[from + t*MaxComms]);
                 help2 = double(Degree[vertex + t*Nodes]);
                 
                 //if (help1 - help2 < 0) Rcout << "help1 - help2  " << help1 - help2 << "  help1  " << help1 << "  help2  " << help2 << std::endl;
                 ratiovalue = ratiovalue - LogFunction(help1 - help2) + LogFunction(help1);
             }
             
             if(DegreeCorrect == 0)
             {
                 help1 = double(CurrentCommStubs[from + t*MaxComms]);
                 help2 = double(Degree[vertex + t*Nodes]);
                 
                 if(CurrentCommVertices[from] > 1) {
                    ratiovalue = ratiovalue - (help1-help2)*log(double(CurrentCommVertices[from]-1));
                    ratiovalue = ratiovalue + help1*log(double(CurrentCommVertices[from]));
                 }
             }
             
             // now we do from/from
             help1 = double(2*CurrentEdgeMatrix[from * MaxComms + from + t*MaxComms*MaxComms]);
             help2 = double(2*SelfEdgeCounter[t] + 2*fromcount);
             
             if (SelfEdgeCounter > 0) {
                 std::cout << "SelfEdgeCounter, t = " << t << ", vertex = " << vertex << ", SEC = " << SelfEdgeCounter[t] << std::endl;
             }
             if (help1 - help2 < 0) std::cout << "help1 - help2 (2)  " << help1 - help2 << std::endl;
             if (help1 < 0) std::cout << "help1 (2) " << help1 << std::endl;
             
             ratiovalue = ratiovalue + .5*LogFunction(help1 - help2) - .5*LogFunction(help1);
             
             // now we add in the terms corresponding to dest
             if(DegreeCorrect == 1)
             {
                 help1 = double(CurrentCommStubs[destination + t*MaxComms]);
                 help2 = double(Degree[vertex + t*Nodes]);
                 ratiovalue = ratiovalue - LogFunction(help1 + help2) + LogFunction(help1);
             }
             
             if(DegreeCorrect == 0)
             {
                 help1 = double(CurrentCommStubs[destination + t*MaxComms]);
                 help2 = double(Degree[vertex + t*Nodes]);
                 ratiovalue = ratiovalue - (help1+help2)*log(double(CurrentCommVertices[destination]+1));
                 if( CurrentCommVertices[destination] > 0)
                    ratiovalue = ratiovalue + help1*log(double(CurrentCommVertices[destination]));
             }
             
             // and now dest/dest

             help1 = double(2*CurrentEdgeMatrix[destination * MaxComms + destination + t*MaxComms*MaxComms]);
             help2 = double(2*SelfEdgeCounter[t] + 2*destcount);
             
             if (help1 + help2 < 0) std::cout << "help1 + help2  (3)" << help1 + help2 << std::endl;
             if (help1 < 0) std::cout << "help1 (3)" << help1 << std::endl;
             
             ratiovalue = ratiovalue + .5*LogFunction(help1 + help2) - .5*LogFunction(help1);
             
             
         }
        
         if (Directed)
         {
             
             //Rcout << "CP_Directed" << std::endl;
             // if the total degree of the vertex is zero we know nothing about it
             // in this case we don't ever change its community
             // at the end we put all degree zeroes into their own group
             if(DegreeCorrect == 1)
             {
                 if(Degree[vertex + t*Nodes] == 0 && outDegree[vertex + t*Nodes] == 0)
                     return 0;
             }
         
             // 1) Communities going into the vertex
             // we first add up all the cross-terms (between communities that are not from / destination)
             for(i=0; i < ActualDiffComms[t]; i++)
             {
                 // we lost NeighborSet[i] to NeighborIndex[i] from the from comm
                 // we gain the same amount in the destination comm
                 // IFF the comms were not from and destination
                 if((NeighborIndex[i + t*MaxComms] != from) && (NeighborIndex[i + t*MaxComms] != destination))
                 {
                     help1 = CurrentEdgeMatrix[NeighborIndex[i + t*MaxComms] * MaxComms + from + t*MaxComms*MaxComms];
                     help2 = CurrentEdgeMatrix[NeighborIndex[i + t*MaxComms] * MaxComms + destination + t*MaxComms*MaxComms];
                     help3 = NeighborSet[i];
             
                     ratiovalue = ratiovalue + LogFunction(help1-help3) - LogFunction(help1);
                     ratiovalue = ratiovalue + LogFunction(help2+help3) - LogFunction(help2);
                 }
             
                 if(NeighborIndex[i + t*MaxComms] == from)
                     fromcount = NeighborSet[i + t*MaxComms];
                     //Rcout << "fromcount_" << fromcount << std::endl;
             
                 if(NeighborIndex[i + t*MaxComms] == destination)
                     destcount = NeighborSet[i + t*MaxComms];
                     //Rcout << "destcount_" << destcount << std::endl;
             }
         
             for(i=0; i < outActualDiffComms[t]; i++)
             {
                 // we lost NeighborSet[i] to NeighborIndex[i] from the from comm
                 // we gain the same amount in the destination comm
                 // IFF the comms were not from and destination
                 if((outNeighborIndex[i + t*MaxComms] != from) && (outNeighborIndex[i + t*MaxComms] != destination))
                 {
                     // do update NOTE: each community mcc' gets updated once if it had edges switch out
                     // which is correct, remembering that mcc' is symmetric (///undirected case) and we only count c < c' here
             
                     help1 = CurrentEdgeMatrix[from * MaxComms + outNeighborIndex[i + t*MaxComms] + t*MaxComms*MaxComms];
                     help2 = CurrentEdgeMatrix[destination * MaxComms + outNeighborIndex[i + t*MaxComms] + t*MaxComms*MaxComms];
                     help3 = outNeighborSet[i + t*MaxComms];
             
                     ratiovalue = ratiovalue + LogFunction(help1-help3) - LogFunction(help1);
                     ratiovalue = ratiovalue + LogFunction(help2+help3) - LogFunction(help2);
                 }
             
                 if(outNeighborIndex[i + t*MaxComms] == from)
                     outfromcount += outNeighborSet[i + t*MaxComms];
             
                 if(outNeighborIndex[i + t*MaxComms] == destination)
                     outdestcount += outNeighborSet[i + t*MaxComms];
             }
        
             // now we add in the term corresponding to from / dest
             help1 = double(CurrentEdgeMatrix[from * MaxComms + destination + t*MaxComms*MaxComms]); /// total white to black
             help2 = double(fromcount-outdestcount); /// white to white - white to black
             help3 = double(CurrentEdgeMatrix[destination * MaxComms + from + t*MaxComms*MaxComms]); /// total black to white
             help4 = double(outfromcount-destcount); /// white to white - black to white
            
             ratiovalue = ratiovalue + LogFunction(help1 + help2) - LogFunction(help1); /// w 2 b - self-edges already excluded
             ratiovalue = ratiovalue + LogFunction(help3 + help4) - LogFunction(help3); /// b 2 w
             
             // now we do from/from
             help1 = double(CurrentEdgeMatrix[from * MaxComms + from + t*MaxComms*MaxComms]);
             help2 = double(SelfEdgeCounter[t] + fromcount + outfromcount);
             
             std::cout << "SelfEdgeCounter, t = " << t << ", SEC = " << SelfEdgeCounter[t] << std::endl;
             if (help1 - help2 < 0) std::cout << "help1 - help2 (2)  " << help1 - help2 << std::endl;
             if (help1 < 0) std::cout << "help1 (2) " << help1 << std::endl;
             
             ratiovalue = ratiovalue + LogFunction(help1 - help2) - LogFunction(help1);
             
             
             // and now dest/dest
             help1 = double(CurrentEdgeMatrix[destination * MaxComms + destination + t*MaxComms*MaxComms]);
             help2 = double(SelfEdgeCounter[t] + destcount + outdestcount);
             ratiovalue = ratiovalue + LogFunction(help1 + help2) - LogFunction(help1);
             
             if(DegreeCorrect == 1)
             {
                 // now we add in the terms corresponding to from
                 ///in
                 help1 = double(CurrentCommEnds[from + t*MaxComms]);
                 help2 = double(Degree[vertex + t*Nodes]);
                 ratiovalue = ratiovalue - LogFunction(help1 - help2) + LogFunction(help1);
                 
                 ///out
                 help1 = double(CurrentCommStubs[from + t*MaxComms]);
                 help2 = double(outDegree[vertex + t*Nodes]);
                 ratiovalue = ratiovalue - LogFunction(help1 - help2) + LogFunction(help1);
                 
                 // now we add in the terms corresponding to dest
                 /// in
                 help1 = double(CurrentCommEnds[destination + t*MaxComms]);
                 help2 = double(Degree[vertex + t*Nodes]);
                 ratiovalue = ratiovalue - LogFunction(help1 + help2) + LogFunction(help1);
                 
                 /// out
                 help1 = double(CurrentCommStubs[destination + t*MaxComms]);
                 help2 = double(outDegree[vertex + t*Nodes]);
                 ratiovalue = ratiovalue - LogFunction(help1 + help2) + LogFunction(help1);
             
             }
         
             if(DegreeCorrect == 0)
             {
                 ///in
                 help1 = double(CurrentCommEnds[from + t*MaxComms]);
                 help2 = double(Degree[vertex + t*Nodes]);
                 if(CurrentCommVertices[from] > 1) {
                     ratiovalue = ratiovalue - (help1 - help2) * log(double(CurrentCommVertices[from]-1)); ///update
                     ratiovalue = ratiovalue + help1 * log(double(CurrentCommVertices[from]));  ///current
                 }
                 
                 ///out
                 help1 = double(CurrentCommStubs[from + t*MaxComms]);
                 help2 = double(outDegree[vertex + t*Nodes]);
                 if(CurrentCommVertices[from] > 1) {
                     ratiovalue = ratiovalue - (help1 - help2) * log(double(CurrentCommVertices[from]-1));
                     ratiovalue = ratiovalue + help1 * log(double(CurrentCommVertices[from]));
                 }
                 
                 ///in
                 help1 = double(CurrentCommEnds[destination + t*MaxComms]);
                 help2 = double(Degree[vertex + t*Nodes]);
                 ratiovalue = ratiovalue - (help1 + help2) * log(double(CurrentCommVertices[destination]+1));
                 ratiovalue = ratiovalue + help1 * log(double(CurrentCommVertices[destination]));
                 
                 ///out
                 help1 = double(CurrentCommStubs[destination + t*MaxComms]);
                 help2 = double(outDegree[vertex + t*Nodes]);
                 ratiovalue = ratiovalue - (help1 + help2) * log(double(CurrentCommVertices[destination]+1));
                 ratiovalue = ratiovalue + help1 * log(double(CurrentCommVertices[destination]));
              }
             
        }
        
    }
        return ratiovalue;
}

void UpdateMatrices(int vertex, int option, int from, int destination)
{
    //Rcout << "Upate Matrices " << std::endl;
    int i;
    double fromcount, destcount, outfromcount, outdestcount;
    
    if ( !Directed )
    {
        
        if(option == 0)
        {
            CurrentCommVertices[from]--;
            CurrentCommVertices[destination]++;
            
            for (int t = 0; t < T; t++) {
                
                fromcount = 0;
                destcount = 0;
                outfromcount = 0;
                outdestcount = 0;
                
                CurrentCommStubs[from + t*MaxComms] -= Degree[vertex + t*Nodes];
                CurrentCommStubs[destination + t*MaxComms] += Degree[vertex + t*Nodes];
            
                for(i=0; i < ActualDiffComms[t]; i++)
                {
                    if((NeighborIndex[i + t*MaxComms] != from) && (NeighborIndex[i + t*MaxComms] != destination))
                    {
                        // do update NOTE: each community mcc' gets updated once if it had edges switch out
                        // which is correct, remembering that mcc' is symmetric and we only count c < c' here
                        CurrentEdgeMatrix[from * MaxComms + NeighborIndex[i + t*MaxComms] + t*MaxComms*MaxComms] -= NeighborSet[i + t*MaxComms];
                        CurrentEdgeMatrix[NeighborIndex[i + t*MaxComms] * MaxComms + from + t*MaxComms*MaxComms] -= NeighborSet[i + t*MaxComms];
                        
                        CurrentEdgeMatrix[destination * MaxComms + NeighborIndex[i + t*MaxComms] + t*MaxComms*MaxComms] += NeighborSet[i + t*MaxComms];
                        CurrentEdgeMatrix[NeighborIndex[i + t*MaxComms] * MaxComms + destination + t*MaxComms*MaxComms] += NeighborSet[i + t*MaxComms];
                    }
                    
                    if(NeighborIndex[i + t*MaxComms] == from)
                        fromcount = NeighborSet[i + t*MaxComms];
                    if(NeighborIndex[i + t*MaxComms] == destination)
                        destcount = NeighborSet[i + t*MaxComms];
                }
                
                CurrentEdgeMatrix[from * MaxComms + from + t*MaxComms*MaxComms] -= (SelfEdgeCounter[t] + fromcount);
                CurrentEdgeMatrix[destination * MaxComms + destination + t*MaxComms*MaxComms] += (SelfEdgeCounter[t] + destcount);
                CurrentEdgeMatrix[from * MaxComms + destination + t*MaxComms*MaxComms] += (fromcount - destcount);
                CurrentEdgeMatrix[destination * MaxComms + from + t*MaxComms*MaxComms] += (fromcount - destcount);
                
            }
        }
        
        if(option == 1)
        {
            BestCommVertices[from]--;
            BestCommVertices[destination]++;
            
            for (int t = 0; t < T; t++) {
                
                fromcount = 0;
                destcount = 0;
                outfromcount = 0;
                outdestcount = 0;
                
                BestCommStubs[from + t*MaxComms] -= Degree[vertex + t*Nodes];
                BestCommStubs[destination + t*MaxComms] += Degree[vertex + t*Nodes];
                
                for(i=0; i < ActualDiffComms[t]; i++)
                {
                    if((NeighborIndex[i + t*MaxComms] != from) && (NeighborIndex[i + t*MaxComms] != destination))
                    {
                        // do update NOTE: each community mcc' gets updated once if it had edges switch out
                        // which is correct, remembering that mcc' is symmetric and we only count c < c' here
                        BestEdgeMatrix[from * MaxComms + NeighborIndex[i + t*MaxComms] + t*MaxComms*MaxComms] -= NeighborSet[i + t*MaxComms];
                        BestEdgeMatrix[NeighborIndex[i + t*MaxComms] * MaxComms + from + t*MaxComms*MaxComms] -= NeighborSet[i + t*MaxComms];
                        
                        BestEdgeMatrix[destination * MaxComms + NeighborIndex[i + t*MaxComms] + t*MaxComms*MaxComms] += NeighborSet[i + t*MaxComms];
                        BestEdgeMatrix[NeighborIndex[i + t*MaxComms] * MaxComms + destination + t*MaxComms*MaxComms] += NeighborSet[i + t*MaxComms];
                    }
                    
                    if(NeighborIndex[i + t*MaxComms] == from)
                        fromcount = NeighborSet[i + t*MaxComms];
                    if(NeighborIndex[i + t*MaxComms] == destination)
                        destcount = NeighborSet[i + t*MaxComms];
                }
                
                BestEdgeMatrix[from * MaxComms + from + t*MaxComms*MaxComms] -= (SelfEdgeCounter[t] + fromcount);
                BestEdgeMatrix[destination * MaxComms + destination + t*MaxComms*MaxComms] += (SelfEdgeCounter[t] + destcount);
                BestEdgeMatrix[from * MaxComms + destination + t*MaxComms*MaxComms] += (fromcount - destcount);
                BestEdgeMatrix[destination * MaxComms + from + t*MaxComms*MaxComms] += (fromcount - destcount);
                
            }
        }
                
        return;
    }
    
    
    if ( Directed )
    {
        
        if(option == 0)
        {
            CurrentCommVertices[from]--;
            CurrentCommVertices[destination]++;
            
            for (int t = 0; t < T; t++)
            {
                fromcount = 0;
                destcount = 0;
                outfromcount = 0;
                outdestcount = 0;
                
                CurrentCommStubs[from + t*MaxComms] -= outDegree[vertex + t*Nodes] ;
                CurrentCommStubs[destination + t*MaxComms] += outDegree[vertex + t*Nodes];
                CurrentCommEnds[from + t*MaxComms] -= Degree[vertex+ t*Nodes];
                CurrentCommEnds[destination + t*MaxComms] += Degree[vertex + t*Nodes];
            
                for(i=0; i < ActualDiffComms[t]; i++)
                {
                    if((NeighborIndex[i] != from) && (NeighborIndex[i] != destination))
                    {
                        CurrentEdgeMatrix[NeighborIndex[i + t*MaxComms] * MaxComms + from + t*MaxComms*MaxComms] -=
                        NeighborSet[i + t*MaxComms];
                        CurrentEdgeMatrix[NeighborIndex[i + t*MaxComms] * MaxComms + destination + t*MaxComms*MaxComms] +=
                        NeighborSet[i + t*MaxComms];
                    }
                    
                    if(NeighborIndex[i + t*MaxComms] == from)
                        fromcount = NeighborSet[i + t*MaxComms];
                    
                    if(NeighborIndex[i + t*MaxComms] == destination)
                        destcount = NeighborSet[i + t*MaxComms];
                    
                }
                
                for(i=0; i < outActualDiffComms[t]; i++)
                {
                    if((outNeighborIndex[i + t*MaxComms] != from) && (outNeighborIndex[i + t*MaxComms] != destination))
                    {
                        CurrentEdgeMatrix[from * MaxComms + outNeighborIndex[i + t*MaxComms] + t*MaxComms*MaxComms] -= outNeighborSet[i + t*MaxComms];
                        CurrentEdgeMatrix[destination * MaxComms + outNeighborIndex[i + t*MaxComms] + t*MaxComms*MaxComms] += outNeighborSet[i + t*MaxComms];
                    }
                    
                    if(outNeighborIndex[i + t*MaxComms] == from)
                        outfromcount = outNeighborSet[i + t*MaxComms];
                    
                    if(outNeighborIndex[i + t*MaxComms] == destination)
                        outdestcount = outNeighborSet[i + t*MaxComms];
                    
                }
                
                CurrentEdgeMatrix[from * MaxComms + from + t*MaxComms*MaxComms] -=
                (fromcount + outfromcount + SelfEdgeCounter[t]);
                CurrentEdgeMatrix[destination * MaxComms + destination + t*MaxComms*MaxComms] +=
                (destcount + outdestcount + SelfEdgeCounter[t]);
                CurrentEdgeMatrix[destination * MaxComms + from + t*MaxComms*MaxComms] +=
                (outfromcount - destcount); /// reminder: self-edges NOT included with NeighborSet
                CurrentEdgeMatrix[from * MaxComms + destination + t*MaxComms*MaxComms] +=
                (fromcount - outdestcount);
            }
        }
        
        if(option == 1)
        {
            BestCommVertices[from]--;
            BestCommVertices[destination]++;
            
            for (int t = 0; t < T; t++)
            {
                fromcount = 0;
                destcount = 0;
                outfromcount = 0;
                outdestcount = 0;
                
                BestCommStubs[from + t*MaxComms] -= outDegree[vertex + t*Nodes] ;
                BestCommStubs[destination + t*MaxComms] += outDegree[vertex + t*Nodes];
                BestCommEnds[from + t*MaxComms] -= Degree[vertex+ t*Nodes];
                BestCommEnds[destination + t*MaxComms] += Degree[vertex + t*Nodes];
                
                
                for(i=0; i < ActualDiffComms[t]; i++)
                {
                    if((NeighborIndex[i + t*MaxComms] != from) && (NeighborIndex[i + t*MaxComms] != destination))
                    {
                        BestEdgeMatrix[NeighborIndex[i + t*MaxComms] * MaxComms + from + t*MaxComms*MaxComms] -= NeighborSet[i + t*MaxComms];
                        BestEdgeMatrix[NeighborIndex[i + t*MaxComms] * MaxComms + destination + t*MaxComms*MaxComms] += NeighborSet[i + t*MaxComms];
                    }
                    
                    if(NeighborIndex[i + t*MaxComms] == from)
                        fromcount = NeighborSet[i + t*MaxComms];
                    
                    if(NeighborIndex[i + t*MaxComms] == destination)
                        destcount = NeighborSet[i + t*MaxComms];
                    
                }
                
                for(i=0; i < outActualDiffComms[t]; i++)
                {
                    if((outNeighborIndex[i + t*MaxComms] != from) && (outNeighborIndex[i + t*MaxComms] != destination))
                    {
                        BestEdgeMatrix[from * MaxComms + outNeighborIndex[i + t*MaxComms] + t*MaxComms*MaxComms] -= outNeighborSet[i + t*MaxComms];
                        BestEdgeMatrix[destination * MaxComms + outNeighborIndex[i + t*MaxComms] + t*MaxComms*MaxComms] += outNeighborSet[i + t*MaxComms];
                    }
                    
                    if(outNeighborIndex[i + t*MaxComms] == from)
                        outfromcount = outNeighborSet[i + t*MaxComms];
                    
                    if(outNeighborIndex[i + t*MaxComms] == destination)
                        outdestcount = outNeighborSet[i + t*MaxComms];
                    
                }
          
                /// already dealt with self edges
                BestEdgeMatrix[from * MaxComms + from + t*MaxComms*MaxComms] -=
                (fromcount + outfromcount + SelfEdgeCounter[t]);
                BestEdgeMatrix[destination * MaxComms + destination + t*MaxComms*MaxComms] +=
                (destcount + outdestcount + SelfEdgeCounter[t]);
                BestEdgeMatrix[destination * MaxComms + from + t*MaxComms*MaxComms] +=
                (outfromcount - destcount); /// reminder: self-edges NOT included with NeighborSet
                BestEdgeMatrix[from * MaxComms + destination + t*MaxComms*MaxComms] +=
                (fromcount - outdestcount);
                
            }
        }
        return;
    }
    
     
}


// This function returns zero if x = 0, otherwise it returns x*log(x)

double LogFunction(double x)
{
    long double tol = 0.0000000001; // this prevents loops due to numerical errors.
    //double tol = std::numeric_limits<double>::epsilon();
    
    if(x < -1.0*tol)
    {
        throw std::range_error("SOMETHING WRONG HAS OCCURRED STOP! x is below zero: ");
    }
    
    if (x > -1.0*tol && x <= tol) {
        //std::cout << "x within tol" << std::endl;
        return 0;
    }
    
    if (x == 0)
        return 0;
    else
        return x*log(x);
    
}

IntegerVector randomComms(const int Nodes, const int MaxComms) {
    RNGScope scope;
    NumericVector randomComms1 = MaxComms * runif(Nodes);
    IntegerVector randomComms2(Nodes);
    int i;
    for (i = 0; i < Nodes; i++)
        randomComms2[i] = (int)std::floor(randomComms1[i]);
    
    return(randomComms2);
    
}


//*********************** KLt FUNCTIONS **************************************************


void Setup(int Nodes, int MaxComms, bool Directed) {  //, const IntegerVector seedComms)
    
    /// Make space in Global Vars
    // TrueState.assign(Nodes, 0); // This records the true communities if they exist read in from the file
    
    Degree.assign(T*Nodes, 0.0);  // Degree of nodes in the network or in Degree if directed
    Count.assign(T*Nodes, 0);  // Degree of nodes in the network or in Degree if directed
    LastEmpty.assign(T*Nodes, 0);
    ActualDiffComms.assign(T, 0);
    
    if (Directed)
    {
        // Directed-only Initializations
        outDegree.assign(T*Nodes, 0.0);
        outCount.assign(T*Nodes, 0);
        outLastEmpty.assign(T*Nodes, 0);
        
        BestCommEnds.assign(T*MaxComms, 0.0); //for directed, keeps tally of edges ending in a class
        // CurrentCommEnds.assign(MaxComms, 0.0); //for directed, keeps tally of edges ending in a class
        
        //Temps don't need to be stored for all time periods simultaneously
        outNeighborIndex.assign(T*MaxComms, 0);
        outTempNeighborIndex.assign(MaxComms, 0);
        outNeighborSet.assign(T*MaxComms, 0.0);
        outTempNeighborSet.assign(MaxComms, 0.0);
        
        outActualDiffComms.assign(T, 0);
    }
    
    //not time-dependent
    BestState.assign(Nodes, 0);
    BestCommVertices.assign(MaxComms, 0); //[MaxComms]
    CurrentState.assign(Nodes, 0);
    CurrentCommVertices.assign(MaxComms, 0); //[MaxComms]
    
    BestCommStubs.assign(T*MaxComms, 0.0); //if directed, keeps tally of edges originating in a class
    BestCommEnds.assign(T*MaxComms, 0.0); //if directed, keeps tally of edges originating in a class
    BestEdgeMatrix.assign(T*MaxComms*MaxComms, 0.0);

    CurrentCommStubs.assign(T*MaxComms, 0.0);
    CurrentCommEnds.assign(T*MaxComms, 0.0);
    CurrentEdgeMatrix.assign(T*MaxComms*MaxComms, 0.0);
    
    //Temps don't need to be stored for all time periods simultaneously
    NeighborIndex.assign(T*MaxComms, 0); // [MaxComms] //lists the comm
    TempNeighborIndex.assign(MaxComms, 0);
    NeighborSet.assign(T*MaxComms, 0.0);  // lists the weight of edges to that comm
    TempNeighborSet.assign(MaxComms, 0.0);
    
}


// [[Rcpp::export]]

List RunKLt(SEXP edgelistTime, const int maxComms, const bool directed, const int klPerNetwork, const bool degreeCorrect, const int nodes)
    {
        
    int i, j, k;
    Rcpp::List EdgelistTime(edgelistTime);
    T = EdgelistTime.size();
    Rcout << T << " Time Stages" << std::endl;
    Nodes = nodes;
    Rcout << Nodes << " Nodes" << std::endl;
    MaxComms = maxComms;
    Directed = directed;
    KLPerNetwork = klPerNetwork;
    DegreeCorrect = degreeCorrect;
    
    Setup(Nodes, MaxComms, Directed);
    
    Rcpp::List AdjListT(T);
    Rcpp::List AdjListWeightT(T);
    Rcpp::List outAdjListT(T);
    Rcpp::List outAdjListWeightT(T);
    
    int t = 0;
    for (List::iterator it = EdgelistTime.begin(); it != EdgelistTime.end(); ++it )
    {
        NumericMatrix edgelist = as<NumericMatrix>(*it);
        int counter = edgelist.nrow();
        
        // Setup Degree
        if (!Directed)
        {
            for(i=0; i < counter; i++) // First we count the degrees by scanning through the list once
            {

                Degree[edgelist(i, 0) + t*Nodes]+=edgelist(i, 2);
                Degree[edgelist(i, 1) + t*Nodes]+=edgelist(i, 2);
                Count[edgelist(i, 0) + t*Nodes]++;
                Count[edgelist(i, 1) + t*Nodes]++;
                
            }
        }
        
        int outmaxCount = 1; //have to set a default for directed.
        
        if (Directed)
        {
            for(i=0; i < counter; i++)
            {
                Degree[edgelist(i, 1) + t*Nodes]+=edgelist(i, 2);
                outDegree[edgelist(i, 0) + t*Nodes]+=edgelist(i, 2);
                Count[edgelist(i, 1) + t*Nodes]++;
                outCount[edgelist(i, 0) +t*Nodes]++;
            }
            
            outmaxCount = *std::max_element(outCount.begin(), outCount.end());
        }
        
        // Setup AdjLists
        // Create big enough adjlist matrices (less efficient than dynamic mem cpp-only implementation, but seems necessary)
        // keep at this scope so it's set correctly
        // undirected & directed
        int maxCount = *std::max_element(Count.begin() + t*Nodes, Count.end());
        IntegerMatrix AdjList(Nodes, maxCount);
        NumericMatrix AdjListWeight(Nodes, maxCount);
        // directed only
        IntegerMatrix outAdjList(Nodes, outmaxCount);
        NumericMatrix outAdjListWeight(Nodes, outmaxCount);

        //Fill in adj matrices
        if (!Directed)
        {
            
            for(i=0; i < counter; i++)
            {
                
                AdjList(edgelist(i, 0), LastEmpty[edgelist(i, 0) + t*Nodes]) = edgelist(i, 1);
                AdjListWeight(edgelist(i, 0), LastEmpty[edgelist(i, 0) + t*Nodes]) = edgelist(i, 2);
                LastEmpty[edgelist(i, 0) + t*Nodes]+=1;
                
                AdjList(edgelist(i, 1), LastEmpty[edgelist(i, 1) + t*Nodes]) = edgelist(i, 0);
                AdjListWeight(edgelist(i, 1), LastEmpty[edgelist(i, 1) + t*Nodes]) = edgelist(i, 2);
                LastEmpty[edgelist(i, 1) + t*Nodes]+=1;
            }
            
            
        }
        
        if (Directed) {
            
            for(i=0; i < counter; i++)
            {
                
                AdjList(edgelist(i, 1), LastEmpty[edgelist(i, 1) + t*Nodes]) = edgelist(i, 0);
                AdjListWeight(edgelist(i, 1), LastEmpty[edgelist(i, 1) + t*Nodes]) = edgelist(i, 2);
                LastEmpty[edgelist(i, 1) + t*Nodes]+=1;
                
                outAdjList(edgelist(i, 0), outLastEmpty[edgelist(i, 0) + t*Nodes]) = edgelist(i, 1);
                outAdjListWeight(edgelist(i, 0), outLastEmpty[edgelist(i, 0) + t*Nodes]) = edgelist(i, 2);
                
                outLastEmpty[edgelist(i, 0) + t*Nodes]+=1;
            }
            
        }
        
        AdjListT[t] = AdjList;
        AdjListWeightT[t] = AdjListWeight;
        
        if (Directed) {
            outAdjListT[t] = outAdjList;
            outAdjListWeightT[t] = outAdjListWeight;
        }
        
        t++;
    }
    
    //////////////////////////////////////// End Setup ////////////////////////////////////////
    
    int MaxIndex;
    double CurrentScore;  // records the current log-likelihood
    int MaxVertex;  // this records the index of the largest vertex ratio found so far
    double MaxRatio;  // records the value of the ratio, actually it's the log of the ratio
    int MaxPriority; // records the community that the vertex wants to go to
    double ProposalRatio, value;
    int Priority;
    std::vector<int> UpdateIndex(Nodes, -1);
    std::vector<int> ChangeSet(Nodes, 0);
        
    // For reporting best state
    std::vector<int> SavedState(Nodes, 0);
    std::vector<int> SavedCommVertices(MaxComms, 0);
    std::vector<double> SavedCommStubs(MaxComms, 0.0);
    std::vector<double> SavedCommEnds(MaxComms, 0.0);
    std::vector<double> SavedEdgeMatrix(MaxComms*MaxComms, 0.0);
    
    long double tolerance = 0.00000001; // this prevents loops due to numerical errors.
    double HighestScore = -std::numeric_limits<double>::max( );
        
    for (int KL = 0; KL < KLPerNetwork; KL++)
    {
        long int tempvertex = 0;
        double prevMaxScore = -std::numeric_limits<double>::max( );
        
        // Builds objects (Vertices, Stubs, Matrix) based on Comm init (currently random init only)
        InitializeKLt(AdjListT, AdjListWeightT);
        
        // Returns log of the initial score:
        MaxScore = ComputeInitialScore();
        
        while(MaxScore >= prevMaxScore + tolerance)
        {
            Rcout << "MAX SCORE IS: " << MaxScore << std::endl;
            // we start with everything equal to the best values
            CurrentScore = MaxScore;
            prevMaxScore = MaxScore;
            MaxIndex = -1;

            // ChangeSet records which vertices are able to move, in that they haven't already moved during this KL step.
            // Update index will tell when the vertex was chosen to move.
            for(i=0; i < Nodes; i++)
            {
                CurrentState[i] = BestState[i]; //BestState set in initialize
                ChangeSet[i] = i; //re-initialized every time
                UpdateIndex[i] = -1;
            }

            for(i=0; i < MaxComms; i++)
            {
                CurrentCommVertices[i] = BestCommVertices[i];
                
                for (t = 0; t < T; t++) {
                    CurrentCommStubs[i + t*MaxComms] = BestCommStubs[i + t*MaxComms];
                    if ( Directed )
                        CurrentCommEnds[i + t*MaxComms] = BestCommEnds[i + t*MaxComms];
                    for(j=0; j < MaxComms; j++)
                        CurrentEdgeMatrix[i*MaxComms+j + t*MaxComms*MaxComms] = BestEdgeMatrix[i*MaxComms+j + t*MaxComms*MaxComms];
                }
            }

            // This loop moves each vertex once
            // Note that we DONT reinitialize changeset as this is unnecessary
            // This would make it a factor of 2 slower.
            for(i=0; i < Nodes; i++)
            {
                MaxVertex = 0;
                MaxRatio = -std::numeric_limits<double>::max( );
                MaxPriority = 0;
                
                // This loop selects which vertex to move;  as i increments up we have one less possible move
                for(j=0; j < Nodes-i; j++)
                {

                    // get proposal and proposal ratio for ChangeSet[j]
                    Priority = 0;
                    ProposalRatio = -std::numeric_limits<double>::max( );
                    // we first compute the neighbor set of the vertex, this is fixed
                    // and the same for every change,
                    // computing this first makes this more efficient
                    // zero indicates run with current communities
                    
                    for (t = 0; t < T; t++) {
                    
                        if (!Directed) {
                            ComputeNeighborSet(ChangeSet[j], 0, AdjListT[t], AdjListWeightT[t], AdjListT[t], AdjListWeightT[t], t);
                        } else {
                            ComputeNeighborSet(ChangeSet[j], 0, AdjListT[t], AdjListWeightT[t], outAdjListT[t], outAdjListWeightT[t], t);
                        }
                    }
                    
                    /*Rcout << " " << std::endl;
                    Rcout << "NeighborSet i,j - "<< i << ", " << j << ": " << std::endl;
                    for (int l = 0; l < T*MaxComms; l++)
                        Rcout << NeighborSet[l] << " _ ";
                    Rcout << std::endl;*/
                        
                    // select how (to which community) to move the vertex
                    for(k=0; k < MaxComms; k++)
                    {
                        // we compute the value of vertex ChangeSet[j] going to k
                        // we DONT allow a vertex to remain where it was
                        // This is essential to enforce so that it will go downhill and not be greedy
                        if(k != CurrentState[ChangeSet[j]])
                        {
                            value = ComputeProposal(ChangeSet[j], CurrentState[ChangeSet[j]], k);
                            if(value > ProposalRatio)
                            {
                                Priority = k;
                                ProposalRatio = value;
                            }
                        }
                    }
                    
                    // check whether its higher than what you already have as the max KL move
                    if(ProposalRatio > MaxRatio)
                    {
                        MaxVertex = j;  // Note this is not the vertex j, but the vertex given by ChangeSet[j]
                        MaxRatio = ProposalRatio;
                        MaxPriority = Priority;
                    }
                }

                // now we move it, first recording the current neighbors so that
                // we can update the matrices properly

                for (t = 0; t < T; t++) {
                
                    if (!Directed) {
                        ComputeNeighborSet(ChangeSet[MaxVertex], 0, AdjListT[t], AdjListWeightT[t], AdjListT[t], AdjListWeightT[t], t);
                    } else {
                        ComputeNeighborSet(ChangeSet[MaxVertex], 0, AdjListT[t], AdjListWeightT[t], outAdjListT[t], outAdjListWeightT[t], t);
                    }
                    
                }
            
                // This updates the matrices to represent the vertices new state
                
                UpdateMatrices(ChangeSet[MaxVertex], 0, CurrentState[ChangeSet[MaxVertex]], MaxPriority);
                
                CurrentState[ChangeSet[MaxVertex]] = MaxPriority;

                // we are using logs so we add the maxratio to the current score for the new score
                CurrentScore = CurrentScore + MaxRatio;

                UpdateIndex[ChangeSet[MaxVertex]] = i;
                // we switch it with the last element of changeset, removing it from further consideration
                // until we have moved the other vertices

                tempvertex = ChangeSet[MaxVertex];
                ChangeSet[MaxVertex] = ChangeSet[Nodes-i-1];
                ChangeSet[Nodes-i-1] = tempvertex;
                
                // now if the new state is better than the previous best state we record this
                // MaxIndex in combination with UpdateIndex
                // telling us where we are in the movement of vertices

                if(CurrentScore > MaxScore)
                {
                    MaxScore = CurrentScore;
                    MaxIndex = i; //tells which move is best
                }
            }
                
            // now we update BestState if a change resulted in a higher maximum
            // by implementing all the changes found above

            // There is a potential for speeding this up here.
            
            if(MaxIndex != -1)
            {
                for(i=0; i < Nodes; i++)
                {
                     // we only make the changes to beststate that happened before or equal to maxindex
                     // no other vertex is updated
                     // fortunately the update order is irrelevant to the final result so
                     // we can just do it this way
                     // Because we force all moves to be different, these updates are all a switch of community
                     if(UpdateIndex[i] <= MaxIndex)
                     {
                         // the option 1 does update corresponding to the best states / matrices
                         for (int t = 0; t < T; t++)
                         {
                             if (!Directed) {
                                 ComputeNeighborSet(i, 1, AdjListT[t], AdjListWeightT[t],
                                                    AdjListT[t], AdjListWeightT[t], t);
                             } else {
                                ComputeNeighborSet(i, 1, AdjListT[t], AdjListWeightT[t],
                                                   outAdjListT[t], outAdjListWeightT[t], t);
                             }
                         }
                             
                         UpdateMatrices(i, 1, BestState[i], CurrentState[i]); // 1 does best matrix update
                         
                         BestState[i] = CurrentState[i];
                     }
                 }
            }
        }
        
        if(MaxScore >= HighestScore)
        {
            HighestScore = MaxScore;
            /*if(TrueCommsAvailable == 1) {
               VIValue = ComputeVI();
               NMIValue = ComputeNMI();
             }*/
            
            for(i=0; i < MaxComms; i++)
                
            {
                SavedCommVertices[i] = BestCommVertices[i];
                SavedCommStubs[i] = BestCommStubs[i];
                if (Directed) SavedCommEnds[i] = BestCommEnds[i];
                for(k=0; k < MaxComms; k++)
                    SavedEdgeMatrix[i*MaxComms+k] = BestEdgeMatrix[i*MaxComms+k];
            }
            for(i=0; i < Nodes; i++)
                SavedState[i] = BestState[i];
        }
    }
    
    for(i=0; i < MaxComms; i++)
    {
        BestCommVertices[i] = SavedCommVertices[i];
        BestCommStubs[i] = SavedCommStubs[i];
        if (Directed) SavedCommEnds[i] = BestCommEnds[i];
        for(k=0; k < MaxComms; k++)
            BestEdgeMatrix[i*MaxComms+k] = SavedEdgeMatrix[i*MaxComms+k];
    }
    
    for(i=0; i < Nodes; i++)
        BestState[i] = SavedState[i];
        
    return List::create(Rcpp::Named("FoundComms") = BestState,//note: correspond to ORDER of IDs
                        Rcpp::Named("EdgeMatrix") = BestEdgeMatrix,
                        Rcpp::Named("HighestScore") = HighestScore);
}

// Calculate Vertices, Stubs, Ends, Matrix based on Comm initiation (currently random initiation only)
void InitializeKLt(List AdjListT, List AdjListWeightT)
{
    
    int i, j;
    int neighbor;
    // double sum;
    IntegerVector randComms = randomComms(Nodes, MaxComms);
    
    
    // initialize - these apply to directed and undirected
    
    for(i=0; i < MaxComms; i++)
        { BestCommVertices[i] = 0;}
    
    for(i=0; i < Nodes; i++)
    {
        BestState[i] = randComms[i];
        /*if(InitializationOption == 1)
         BestState[i] = TrueState[i];*/
        BestCommVertices[BestState[i]]++;
    }
    
    
    for (int t = 0; t < T; t ++)
    {
        
        if ( !Directed )
        {
            
            // initialize
            for(i=0; i < MaxComms; i++)
            {
                BestCommStubs[i+t*MaxComms] = 0;
                for(j=0; j < MaxComms; j++)
                    BestEdgeMatrix[i * MaxComms + j + t*MaxComms*MaxComms] = 0;
            }
            
            for(i=0; i < Nodes; i++)
            {
                BestCommStubs[BestState[i] + t*MaxComms] += Degree[i + t*Nodes]; ///stubs will add to twice the total edge weight
            }
            
            /// sum = 0;
            
            for(i=0; i < Nodes; i++)
            {
                for(j=0; j < LastEmpty[i + t*Nodes]; j++)
                {
                    NumericMatrix tmp = AdjListT[t];
                    neighbor = tmp(i, j); ///each edge listed twice, once for each end
                    NumericMatrix tmpw = AdjListWeightT[t];
                    double neighborWeight = tmpw(i, j);
                    /// sum += neighborWeight;
                    
                    /// twice the diagonal of BestEdgeMatrix + everything else in the matrix once should equal = 2x edges
                    
                    if (BestState[neighbor] == BestState[i])
                        neighborWeight = (neighborWeight/2.0);
                    
                    BestEdgeMatrix[(BestState[i])*MaxComms + BestState[neighbor] + t*MaxComms*MaxComms] += neighborWeight;
                    
                }
            }
        }
        
        if ( Directed )
        {
            
            for(i=0; i < MaxComms; i++) /// initialize
            {
                BestCommStubs[i + t*MaxComms] = 0;
                BestCommEnds[i + t*MaxComms] = 0;
                for(j=0; j < MaxComms; j++)
                    BestEdgeMatrix[i*MaxComms+j + t*MaxComms*MaxComms] = 0;
            }
            
            
            for(i=0; i < Nodes; i++)
            {
                
                BestCommStubs[BestState[i] + t*MaxComms] += outDegree[i +  t*Nodes];
                BestCommEnds[BestState[i] + t*MaxComms] += Degree[i + t*Nodes];
            }
            
            // sum = 0;
            
            for(i=0; i < Nodes; i++)
            {
                for(j=0; j < LastEmpty[i + t*Nodes]; j++)
                {
                    NumericMatrix tmp = AdjListT[t];
                    neighbor = tmp(i, j); ///each edge listed twice, once for each end
                    NumericMatrix tmpw = AdjListWeightT[t];
                    
                    BestEdgeMatrix[BestState[neighbor] * MaxComms + BestState[i] + t*MaxComms*MaxComms] += tmpw(i, j);
                    // sum += AdjListWeight(i, j); /// Note self-edges get counted ONCE because we're only interating over out-edges
                    
                    /// count self edges twice?
                    if (i == neighbor && SelfTwice == 2) {
                        int position = BestState[i]*MaxComms + BestState[neighbor];
                        BestEdgeMatrix[position + t*MaxComms*MaxComms] += tmpw(i, j);
                        // sum += AdjListWeight(i, j);
                    }
                    
                    
                    
                }
                
            }
        }
        
    }
    /// also, need to double-check that sum is correct so everything adds to one
    
}
