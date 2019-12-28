int import_lattice(vector<vector<double>> rules, maps<int,vector<vector<double>>> nns){
    ifstream latfile;
    float scale;
    string name, nextline;
    bool again, blank, fileopen = false;
    vector<double> rule, nn;
    int lcount = 0, pointiter, nncount = 0;

    //Prompt for file name
    do{
        printf("Enter the lattice file name: ");
        getline(cin, nextline);
        name = split(nextline, ' ')[0];

        if(!name.empty()) latfile.open(name);

        if(!latfile.is_open()){
            again = yesno("No file was read. Try again? ");
            if(!again) return 0;
        }
    }while(!latfile.is_open());

    //Process rules until a blank line is reached
    blank = false;
    while(!latfile.eof()){
        lcount ++;
        latfile.getline(nextline);
        if(nextline.empty()) break;

        rule = getdoubles(split(nextline, ' '));
        if(rule.size() < 3){
            fprintf(stderr, "Too few numbers in rule on line line %d\n",lcount);
            latfile.close();
            return 0;
        }
        rules.push_back(rule);
    }
    
    //Read nearest neighbor rules
    for(vector<double> nextrule : rules){
       ruleiter ++;
       for(pointiter = 1; pointiter <= nextrule.size() - 2; pointiter ++){
           vector<vector<double>> nextset;
           while(!latfile.eof()){
               lcount ++;
               getline(latfile, nextline);
               if(nextline.empty()) break;
               nn = getdoubles(split(nextline, ' '));
               if(nn.size() != 2) fprintf(stderr, "Insufficient information for nearest neighbor rule on line %d.\n", lcount);
               else{
                   nextset.push_back(nn);
               }
           }
           nns.insert(pair<int, vector<vector<double>>>(nncount++, nextset));
       }
    }


    latfile.close();
    return 1;
}
