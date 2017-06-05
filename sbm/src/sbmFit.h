

IntegerVector randomComms(const int Nodes, const int MaxComms) {
    RNGScope scope;
    NumericVector randomComms1 = MaxComms * runif(Nodes);
    IntegerVector randomComms2(Nodes);
    int i;
    for (i = 0; i < Nodes; i++)
        randomComms2[i] = (int)std::floor(randomComms1[i]);
    
    return(randomComms2);
    
}

