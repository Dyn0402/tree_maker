
void readMuDst(TString InputFileList, TString OutputDir, int energy)
{
    
    Int_t nFiles = 1e+4;
    
    // Load libraries
    gROOT ->Macro("loadMuDst.C");
    gSystem->Load("St_base");
    gSystem->Load("StChain");
    gSystem->Load("StUtilities");
    gSystem->Load("StIOMaker");
    gSystem->Load("StarClassLibrary");
    gSystem->Load("StEvent");
    gSystem->Load("MyAnalysisMaker") ;
    //    gSystem->Load("StRefMultCorr");
    gSystem->Load("StBTofUtil");
    
    // List of member links in the chain
    StChain*                    chain  =  new StChain ;
    StMuDstMaker*          muDstMaker  =  new StMuDstMaker(0,0,"",InputFileList,"MuDst",nFiles) ;
    MyAnalysisMaker*    AnalysisCode   =  new MyAnalysisMaker(muDstMaker) ;

    // Turn off everything but Primary tracks in order to speed up the analysis and eliminate IO
    muDstMaker -> SetStatus("*",0) ;                // Turn off all branches
    muDstMaker -> SetStatus("MuEvent",1) ;          // Turn on the Event data (esp. Event number)
    muDstMaker -> SetStatus("PrimaryTracks",1) ;    // Turn on the primary track data
    muDstMaker -> SetStatus("BTofHit",1) ;
    muDstMaker -> SetStatus("BTofHeader",1) ;
    
    AnalysisCode -> SetOutputFileName(OutputDir) ;
    AnalysisCode -> SetEnergy(energy) ;
    
    Int_t nEvents = 1e5;
    
    nEvents = muDstMaker->chain()->GetEntries();
    
    cout<<"\n############################ Total Event in chain = "<<nEvents<<" ############################\n "<<endl;
    
    // Loop over the links in the chain
    
    chain -> Init() ;

    int iret;
    
    for(int i=0;i<nEvents;i++)
    {

    	cout << "Pre chain->Make(" << i << ")" << endl;
        iret = chain->Make(i);
        cout << "Post chain->Make(" << i << ")" << endl;
        
        if(i%((int)nEvents/50)==0)
            cout<<" => Processing event# = "<<i<<" with return code = "<<iret<<endl;
        
        if(iret){
            cout<<"....File Ends, Analysis DONE...."<<endl; break;
        }
        
    }
    
    chain -> Finish() ;
    
    
    // Cleanup
    delete chain ;
    
}

