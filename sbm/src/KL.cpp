// Called by runKL.R

//*********************** GLOBAL VARIABLES *******************************************************************

//std::vector<int> NeighborIndex; // [MaxComms] //lists the comm
//std::vector<int> outNeighborIndex;
//std::vector<int> TempNeighborIndex;
//std::vector<int> outTempNeighborIndex;

//std::vector<double> NeighborSet;  // lists the weight of edges to that comm
//std::vector<double> outNeighborSet;
//std::vector<double> TempNeighborSet;
//std::vector<double> outTempNeighborSet;

std::vector<double> Degree;  // Degree of nodes in the network or in Degree if directed
std::vector<double> outDegree;

std::vector<int> BestState;
std::vector<int> BestCommVertices; //[MaxComms]
std::vector<double> BestCommStubs; //if directed, keeps tally of edges originating in a class
std::vector<double> BestCommEnds; //for directed, keeps tally of edges ending in a class
std::vector<double> BestEdgeMatrix;

double MaxScore = 0;
bool Directed;
double SelfEdgeCounter; // records the *weight* of self-edges (not doubled) so they can be counted correctly.



//*********************** FUNCTION DECLARATIONS **************************************************************

# functions called by RunKL
# Initialize -- calls Degree
# ComputeInitialScore
# ComputeNeighborSet
# UpdateMatrices
# ComputeProposal

void Setup();
void Initialize(IntegerMatrix AdjList, NumericMatrix AdjListWeight); // initializes the data structures for KL
void RunKL( __ fix __ ); // runs Kernighan-Lin once.
#double ComputeInitialScore(); // computes the initial score after initialization
#void ComputeNeighborSet(int vertex, int option, IntegerMatrix AdjList, NumericMatrix AdjListWeight, IntegerMatrix outAdjList, NumericMatrix outAdjListWeight);  // computes the neighbor set for a vertex, using either best or currentstates
#double ComputeProposal(int vertex, int from, int destination); // computes the value of the particular change
#void UpdateMatrices(int vertex, int option, int from, int destination); // this updates either the best
#double LogFunction(double x); // this returns x*log(x) and zero if x=0
#IntegerVector randomComms(int Nodes, int MaxComms);


void Setup() {  //, const IntegerVector seedComms)
    
    /// Make space in Global Vars
    // ChangeSet.assign(Nodes, 0);
    // UpdateIndex.assign(Nodes, 0);
    // TrueState.assign(Nodes, 0); // This records the true communities if they exist read in from the file
    
    //Degree.assign(Nodes, 0.0);  // Degree of nodes in the network or in Degree if directed
    //Count.assign(Nodes, 0);  // Degree of nodes in the network or in Degree if directed
    //LastEmpty.assign(Nodes, 0);
    
    //Is this necessary?
    BestState.assign(Nodes, 0);
    BestCommVertices.assign(MaxComms, 0); //[MaxComms]
    BestCommStubs.assign(MaxComms, 0.0); //if directed, keeps tally of edges originating in a class
    BestCommEnds.assign(MaxComms, 0.0); //if directed, keeps tally of edges originating in a class
    BestEdgeMatrix.assign(MaxComms*MaxComms, 0.0);
    
    //NeighborIndex.assign(MaxComms, 0); // [MaxComms] //lists the comm
    //TempNeighborIndex.assign(MaxComms, 0);
    //NeighborSet.assign(MaxComms, 0.0);  // lists the weight of edges to that comm
    //TempNeighborSet.assign(MaxComms, 0.0);
    
}

void KL(IntegerMatrix AdjList, NumericMatrix AdjListWeight, IntegerMatrix outAdjList, NumericMatrix outAdjListWeight,
        const bool directed, const int nodes, const int maxComms,
        const NumericVector degree, const NumericVector outdegree)
{
    
    Nodes = nodes;
    Directed = directed;
    MaxComms = maxComms;
    Degree = degree;
    outDegree = outdegree;
    
    Setup();
    
    int i,j,k;
    int MaxIndex = 1;
#double CurrentScore;  // records the current log-likelihood
    int MaxVertex;  // this records the index of the largest vertex ratio found so far
    double MaxRatio;  // records the value of the ratio, actually it's the log of the ratio
    int MaxPriority; // records the community that the vertex wants to go to
    long int tempvertex = 0;
    
    std::vector<int> UpdateIndex;
    std::vector<int> ChangeSet;
    
    std::vector<int> CurrentState;
    std::vector<int> CurrentCommVertices; //[MaxComms]
    std::vector<double> CurrentCommStubs;
    std::vector<double> CurrentCommEnds;
    std::vector<double> CurrentEdgeMatrix;
    
    double prevMaxScore = -std::numeric_limits<double>::max( );
    long double tolerance = 0.00000001; // this prevents loops due to numerical errors.
    
    double ProposalRatio;
    double value;
    int Priority;
    
    Initialize(AdjList, AdjListWeight); // Calculate Vertices, Stubs, Ends, Matrix based on Comm initiation (currently random initiation only)
    
    // This returns the log of the initial score
    /*MaxScore = ComputeInitialScore();
     
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
     CurrentState[i] = BestState[i]; //BestState set in initialize
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
     ComputeNeighborSet(ChangeSet[j], 0, AdjList, AdjListWeight, outAdjList, outAdjListWeight);
     
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
     
     ComputeNeighborSet(ChangeSet[MaxVertex], 0, AdjList, AdjListWeight, outAdjList, outAdjListWeight);
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
     ComputeNeighborSet(i, 1, AdjList, AdjListWeight, outAdjList, outAdjListWeight);
     UpdateMatrices(i, 1, BestState[i], CurrentState[i]); // 1 does best matrix update
     
     BestState[i] = CurrentState[i];
     }
     }
     }
     }*/
    
    return;
}

// Calculate Vertices, Stubs, Ends, Matrix based on Comm initiation (currently random initiation only)
void InitializeKL(IntegerMatrix AdjList, NumericMatrix AdjListWeight)
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
    
    
}

