void test(){
	
	float num[5]={2.02,40.77,2,6.10,8.309};
  
	TH1F* h= new TH1F("","",1000,1,55);
	int size=5;
	int find;
	float get;

	for(int i=0; i<size; i++){
			

		h->Fill(num[i]);
		find = h->FindBin(num[i]);
		cout<<"getbin "<<find<<" num["<<i<<"] "<<num[i]<<endl;
		get=h->GetBinContent(find);
		cout<<"get "<<get<<endl;


	
	}
		
			for(int i=0; i<size; i++){
	
	
	}
	
}
