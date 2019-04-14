// TDD-SBM

#include <Rcpp.h>
#include "sbm.h"
using namespace Rcpp;

// ************************************************************************************************************************

// Prototype functions for time-dependent models
void Setup(int Nodes, int MaxComms, bool Directed);
void InitializeKLt(List AdjList, List AdjListWeight); // initializes the data structures for KL
List sbmtFit(SEXP edgelistTime, const IntegerMatrix & edgelistTotal, const NumericVector & weightsTotal, const int maxComms, const bool directed, const int klPerNetwork, const int degreeCorrect, const int nodes);

// ************************************************************************************************************************

// [[Rcpp::export]]

List sbmtFit(SEXP edgelistTime, const int maxComms, const bool directed, const int klPerNetwork, const int degreeCorrect, const int nodes, const long double tolerance)
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
    Tolerance = tolerance;
    
    Setup(Nodes, MaxComms, Directed);
    
    Rcpp::List AdjListT(T);
    Rcpp::List AdjListWeightT(T);
    Rcpp::List outAdjListT(T);
    Rcpp::List outAdjListWeightT(T);
    
    int t = 0;
    int outmaxCount;
    NumericMatrix edgelist;
    int counter;
    for (List::iterator it = EdgelistTime.begin(); it != EdgelistTime.end(); ++it )
    {
        edgelist = as<NumericMatrix>(*it);
        counter = edgelist.nrow();
        
        // Setup Degree
        if (!Directed)
        {
            for(i=0; i < counter; i++) // First we count the degrees by scanning through the list once
            {
                
                Degree[edgelist(i, 0) + t*Nodes]+=edgelist(i, 2);
                Degree[edgelist(i, 1) + t*Nodes]+=edgelist(i, 2);
                Count[edgelist(i, 0) + t*Nodes]++;
                Count[edgelist(i, 1) + t*Nodes]++;
                if (edgelist(i, 0) == edgelist(i, 1)) {
                    SelfEdgeCounter[t*Nodes + edgelist(i, 0)] += edgelist(i, 2);
                    
                }
            }
        }
        
        outmaxCount = 1; //have to set a default for directed.
        
        if (Directed)
        {
            for(i=0; i < counter; i++)
            {
            
                Degree[edgelist(i, 1) + t*Nodes]+=edgelist(i, 2);
                outDegree[edgelist(i, 0) + t*Nodes]+=edgelist(i, 2);
                Count[edgelist(i, 1) + t*Nodes]++;
                outCount[edgelist(i, 0) +t*Nodes]++;
                if (edgelist(i, 0) == edgelist(i, 1)) {
                    SelfEdgeCounter[t*Nodes + edgelist(i, 0)] += edgelist(i, 2);
                }
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

    Rcout << "End Setup of sbmtFit\n" << std::endl;
    //////////////////////////// End Setup of sbmtFit /////////////////////////////////

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
    double HighestScore = -std::numeric_limits<double>::max( );

    for (int KL = 0; KL < KLPerNetwork; KL++)

    {
        Rcout << "KL " << KL << std::endl;
        long int tempvertex = 0;
        double prevMaxScore = -std::numeric_limits<double>::max( );
        
        // Builds objects (Vertices, Stubs, Matrix) based on Comm init (currently random init only)
        InitializeKLt(AdjListT, AdjListWeightT);
        
        // Returns log of the initial score:
        MaxScore = ComputeInitialScore();
        Rcout << "InitialScore: " << MaxScore << std:: endl;
        
        while(MaxScore >= prevMaxScore + Tolerance)
        {
            Rcout << "MAX SCORE IS: " << MaxScore << std::endl;
            
            /*Rcout << "CurrentCommStubsTotal: " << std::endl;
             for (i = 0; i < MaxComms; i++)
             Rcout << CurrentCommStubsTotal[i] << " " <<std::endl;
             Rcout << "BestCommStubsTotal: " << std::endl;
             for (i = 0; i < MaxComms; i++)
             Rcout << BestCommStubsTotal[i] << " " << std::endl;*/
            
            
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
                
                if (DegreeCorrect == 2 | DegreeCorrect == 3) {
                    CurrentCommStubsTotal[i] = BestCommStubsTotal[i];
                    if (Directed){
                        CurrentCommEndsTotal[i] = BestCommEndsTotal[i];
                    }
                }
                
                
                for (t = 0; t < T; t++) {
                    CurrentCommStubs[i + t*MaxComms] = BestCommStubs[i + t*MaxComms];
                    if ( Directed ) {
                        CurrentCommEnds[i + t*MaxComms] = BestCommEnds[i + t*MaxComms];
                    }
                    for (j=0; j < MaxComms; j++)
                        CurrentEdgeMatrix[i*MaxComms+j + t*MaxComms*MaxComms] = BestEdgeMatrix[i*MaxComms+j + t*MaxComms*MaxComms];
                }
            }
            
            // This loop moves each vertex once
            // Note that we DONT reinitialize changeset as this is unnecessary
            // This would make it a factor of 2 slower.
            for(i=0; i < Nodes; i++)
            {
                //Rcout << "Node i " << i << " ";
                MaxVertex = 0;
                MaxRatio = -std::numeric_limits<double>::max( );
                MaxPriority = 0;
                
                // This loop selects which vertex to move;  as i increments up we have one less possible move
                for(j=0; j < Nodes-i; j++)
                {
                    //Rcout << ", j " << j << " ";
                    //Rcout << "Node j " << j << std::endl;
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
                    // Rcout << "CurrentScore - MaxScore " << CurrentScore - MaxScore <<std::endl;
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
        
        for(i=0; i < MaxComms; i++)
        {
            BestCommVertices[i] = SavedCommVertices[i];
            BestCommStubs[i] = SavedCommStubs[i];
            if (Directed) SavedCommEnds[i] = BestCommEnds[i];
            for(k=0; k < MaxComms; k++)
                BestEdgeMatrix[i*MaxComms+k] = SavedEdgeMatrix[i*MaxComms+k];
        }
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
    
    for(i=0; i < MaxComms; i++) {
        BestCommVertices[i] = 0;
        BestCommStubsTotal[i] = 0;
        if (Directed)
            BestCommEndsTotal[i] = 0;
    }
    
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
                BestCommStubsTotal[BestState[i]] += Degree[i + t*Nodes];
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
                
                BestCommStubsTotal[BestState[i]] += outDegree[i + t*Nodes];
                BestCommEndsTotal[BestState[i]] += Degree[i + t*Nodes];
            }
            
            // sum = 0;
            
            for(i=0; i < Nodes; i++)
            {
                for(j=0; j < LastEmpty[i + t*Nodes]; j++)
                {
                    NumericMatrix tmp = AdjListT[t];
                    neighbor = tmp(i, j); ///each edge listed twice, once for each end
                    NumericMatrix tmpw = AdjListWeightT[t];
                    
                    BestEdgeMatrix[(BestState[neighbor]) * MaxComms + BestState[i] + t*MaxComms*MaxComms] += tmpw(i, j);
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

void Setup(int Nodes, int MaxComms, bool Directed)

{  //, const IntegerVector seedComms)
    
    // Make space in Global Vars
    // TrueState.assign(Nodes, 0); // This records the true communities if they exist read in from the file
    
    SelfEdgeCounter.assign(T*Nodes, 0.0);
    
    Degree.assign(T*Nodes, 0.0);  // Degree of nodes in the network or in Degree if directed
    Count.assign(T*Nodes, 0);  // Degree of nodes in the network or in Degree if directed
    LastEmpty.assign(T*Nodes, 0);
    
    BestCommStubs.assign(T*MaxComms, 0.0); //if directed, keeps tally of edges originating in a class
    BestEdgeMatrix.assign(T*MaxComms*MaxComms, 0.0);
    
    CurrentCommStubs.assign(T*MaxComms, 0.0);
    CurrentEdgeMatrix.assign(T*MaxComms*MaxComms, 0.0);
    
    NeighborSet.assign(T*MaxComms, 0.0);  // lists the weight of edges to that comm
    
    if (Directed)
    {
        // Directed-only Initializations
        outDegree.assign(T*Nodes, 0.0);
        outCount.assign(T*Nodes, 0);
        outLastEmpty.assign(T*Nodes, 0);
        
        CurrentCommEnds.assign(T*MaxComms, 0.0);
        BestCommEnds.assign(T*MaxComms, 0.0); //for directed, keeps tally of edges ending in a class
        
        outNeighborSet.assign(T*MaxComms, 0.0);
        
        CurrentCommEndsTotal.assign(MaxComms, 0.0); //for time-independent degree correction on time-dependent data
        BestCommEndsTotal.assign(MaxComms, 0.0);
        
    }
    
    //not time-dependent
    BestState.assign(Nodes, 0);
    BestCommVertices.assign(MaxComms, 0); //[MaxComms]
    CurrentState.assign(Nodes, 0);
    CurrentCommVertices.assign(MaxComms, 0); //[MaxComms]
    
    CurrentCommStubsTotal.assign(MaxComms, 0.0);
    BestCommStubsTotal.assign(MaxComms, 0.0);
    
}
