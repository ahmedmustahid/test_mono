
#include <vector>
#include <chrono>

using namespace chrono;

    
void weight(float *a, float *ctheta, float *tol, float *wt,int size, float *b);
void fill_hist(float *a, float *ctheta, float *tol, TH2F* h, TH2F* h1);
//void fill_hist(float *a,float *ctheta, float *tol, TH2F* h, TH2F* h1);
void init_ar(float *a, int size);
void create_edge( float *arr, float *tol, float *edges_a);
float dsig(float e, float cthta, float tol1);
int nung_count(float ener, float cangle, float *a, float *b,int size, float tol1);
void tol_e(int nung_num,float e, float cthta,float *a, float *b,int size,  float *tol);
float tol_least(float *tol, int tol_size);

int nung_tot=0;

void efficiency_new_TGraph(){
	auto t1 = high_resolution_clock::now();

		
        TFile* fin = new TFile("l5_500GeV.nung_1.root");
	
        TTree* data = static_cast<TTree*>(fin->Get("evtdata"));
	
	  	
//	TH2D* htheta_sigma = new TH2D("htheta_sigma",";#theta;#sigma of fit",1000,0,3,100,0.2,0.3); 	
	TH1D* hmcp_pfo_theta = new TH1D("hmcp_pfo_theta","Signal Theta Distribution;#theta_{mcp}",100,-1,4); 	
    
	TH1D* hresol_e = new TH1D("hresol_e",";E_{pfo}-E_{mcp}",20,-6,6); 	
	TH1D* hdsig = new TH1D("hdsig",";diff cross (nb)",100,0,2); 	
	TH1D* hresol_phi = new TH1D("hresol_phi",";#phi_{pfo}-#phi_{mcp}",20,-4,4); 	
	TH1D* hresol_theta = new TH1D("hresol_theta",";#theta_{pfo}-#theta_{mcp}",20,-1,1); 	
	TH1D* hmcp_pfo_phi = new TH1D("hmcp_pfo_phi",";#phi_{pfo}-#phi_{mcp}",100,-7,7); 	

    TH1D* hmcp_theta = new TH1D("hmcp_theta",";#theta_{mcp}",100,-1,4); 	
    TH2D* hmcp_e_cos = new TH2D("hmcp_e_cos","nung;cos#theta_{mcp};E_{mcp}",100,-1,1,100,0,250); 
    TH1D* hmcp_pfo_e = new TH1D("hmcp_pfo_e",";E_{mcp}",1000,0,300); 	
//        TH1D* hpfo_theta = new TH1D("hpfo_theta",";#theta;",100,0,4);  
//        TH1D* hpfo_phi = new TH1D("hpfo_phi",";#phi;",100,-3,3);  
//	TCanvas* cmcp_pfo_e = new TCanvas("example","",600,400);
//	TCanvas* ctheta_sigma = new TCanvas("example","",600,400);
	//TCanvas* cpfo_theta_phi = new TCanvas("phi theta","",900,600);
	//cpfo_theta_phi->Divide(2,1);
		
	TH1F* he_nung= new TH1F("he_nung",";E",1000,0,250);
 	TF1* fmcp_pfo_e = new TF1("fmcp_pfo_e","gaus(0)",-50,50);

    	const int NMAX = 10000;
    	int npfos, nmcps, mcp_pdg[NMAX], mcp_index[NMAX], mcp_nparents[NMAX], mcp_parentIndex[NMAX][10], mcp_daughterIndex[NMAX][10],mcp_ndaughters[NMAX], pfo_pdg[NMAX];
    	int j0=0,s=500*500,spn=0,mx=120;
    	float alpha=1.0/137.0, ke=0.6, sig_an=7.0*pow(10,-3);//unit of cross sect in nanobarn
    	float dsig_const = (alpha * ke * sig_an * s * pow(4,j0)*pow((2*spn+1),2))/(64.0 * 3.1416);	
    	cout<<"dsig_const "<<dsig_const<<endl;
    	cout<<"alpha "<<alpha<<endl;
    	float pfo_e[NMAX], mcp_e[NMAX], mcp_theta[NMAX], sigma[NMAX], pfo_phi[NMAX], pfo_theta[NMAX], mcp_phi[NMAX], cangle[NMAX], temp = 1,dm_cs=0.0,dm_cs1=0.0;

 
		

	bool bPassed;

    data->SetBranchAddress("mcp_pdg",mcp_pdg);
    data->SetBranchAddress("pfo_pdg",pfo_pdg);
    data->SetBranchAddress("mcp_e",mcp_e);
  	data->SetBranchAddress("mcp_index",mcp_index);
  	data->SetBranchAddress("mcp_nparents",mcp_nparents);
  	data->SetBranchAddress("mcp_parentIndex",mcp_parentIndex);
  	data->SetBranchAddress("mcp_ndaughters",mcp_ndaughters);
  	data->SetBranchAddress("mcp_daughterIndex",mcp_daughterIndex);
    data->SetBranchAddress("pfo_e",pfo_e);
  	data->SetBranchAddress("pfo_theta",pfo_theta);
  	data->SetBranchAddress("pfo_phi",pfo_phi);
  	data->SetBranchAddress("mcp_theta",mcp_theta);
  	data->SetBranchAddress("mcp_phi",mcp_phi);
  	//data->SetBranchAddress("mcp_theta",mcp_theta);
//  	data->SetBranchAddress("mcp_endx",mcp_endx);
//  	data->SetBranchAddress("mcp_endy",mcp_endy);
//  	data->SetBranchAddress("mcp_endz",mcp_endz);
	data->SetBranchAddress("npfos",&npfos);
	data->SetBranchAddress("nmcps",&nmcps);

	


	int nevents = data->GetEntries();
	nevents = 500;
    	float a[nevents],b[nevents],c[nevents];//for saving e, theta, phi of mcp
	
	//initializing arrays
	for(int i=0; i<nevents; i++){
	
		a[i]= 1000;
		b[i]= 1000;
		c[i]= 1000;
	
	}


    
    for(int ev=0; ev<nevents; ev++){
		
        data->GetEntry(ev);

	//mcp signal
	
	for(int ip=0; ip < nmcps; ip++){
       	
       		//hmcp_e_cos_dsig->Fill(mcp_e[ip],TMath::Cos(mcp_theta[ip]));
		//int bin_no= hmcp_e_cos_dsig->FindBin(mcp_e[ip],TMath::Cos(mcp_theta[ip]));
		//int bin_cont= hmcp_e_cos_dsig->GetBinContent(bin_no);

		//cout<<"bin_no "<<bin_no<<endl;
		//cout<<"bin_cont "<<bin_cont<<endl;
            int m = 0; //m makes sure parent ID are electron positron

            if(mcp_pdg[ip]==22 &&  mcp_nparents[ip]==2){
            
                 int j=0;
                if( ( (mcp_pdg[mcp_parentIndex[ip][j]]==11 &&  mcp_pdg[mcp_parentIndex[ip][j+1]]==-11) ||  (mcp_pdg[mcp_parentIndex[ip][j]]==-11 &&  mcp_pdg[mcp_parentIndex[ip][j+1]]==11) ) ) m++;

            
                int g=0, v=0, vbar=0;

                if(m==1 && mcp_ndaughters[mcp_parentIndex[ip][0]]==3){
            
                    
                    for(int k=0; k< mcp_ndaughters[mcp_parentIndex[ip][0]]; k++){
                        
                         if( mcp_pdg[mcp_daughterIndex[mcp_parentIndex[ip][0]][k]]==22 ) g++;
                         if( mcp_pdg[mcp_daughterIndex[mcp_parentIndex[ip][0]][k]]==12 ) v++;
                         if( mcp_pdg[mcp_daughterIndex[mcp_parentIndex[ip][0]][k]]==-12 ) vbar++;
                                
                    }
                         
                     if(g==1 && v== 1 && vbar==1){ 
                        hmcp_theta->Fill(mcp_theta[ip]);
			

                        
                        a[nung_tot]=mcp_e[ip];
                        b[nung_tot]=mcp_theta[ip];
                        c[nung_tot]=mcp_phi[ip];
                        hmcp_pfo_e->Fill(mcp_e[ip]);
                        nung_tot++;


			hmcp_e_cos->Fill(TMath::Cos(mcp_theta[ip]),mcp_e[ip]);
			//int bin_no=hmcp_e_cos->FindBin(TMath::Cos(mcp_theta[ip]),mcp_e[ip]);
			//cout<<"bin_no "<<bin_no<<endl;
                        if(ip!=8) cout<<"mcp index "<<ip<<endl;
                    }
                }
            }
        }
    	


	

	//pfo_resolution	
    
	int a_size=sizeof(a)/sizeof(a[0]);
        int i=0;
        while(i<a_size && i<nung_tot){ //a is e, b is theta, c is phi->1000 determines the last element
                    
            
            for(int ip=0; ip< npfos; ip++){
                bool e = pfo_e[ip]<=(a[i]+5) && pfo_e[ip]>=(a[i]-5) ;//5 has been determined from mcr-pfo resol plots

                bool theta =  pfo_theta[ip]<=(b[i]+0.2) && pfo_theta[ip]>=(b[i]-0.2);
                bool phi= pfo_phi[ip]<=(c[i]+2) && pfo_phi[ip]>=(c[i]-2);
                bool gamma = pfo_pdg[ip]==22;
		//add some more condiotns

                //int count=0;

                if(e && theta && phi && gamma){

                    
                    hmcp_pfo_theta->Fill(b[i]);
                    hresol_e->Fill(pfo_e[ip]-a[i]);
                    hresol_theta->Fill(pfo_theta[ip]-b[i]);
                    hresol_phi->Fill(pfo_phi[ip]-c[i]);
                    break;// this is not required. change later
                }
            }
        i++;
        }
    }
    	
    	TCanvas* test1= new TCanvas("test1","",900,900);

	test1->cd();
	hmcp_e_cos->Draw();

		
	
	
	cout<<" nung_tot "<<nung_tot<<endl;	
	int a_size=sizeof(a)/sizeof(a[0]);
	float tol[nung_tot];	
	for(int i=0; i<nung_tot;i++){
	
		tol[i]=0;
	}
	int i=0;
	while(i<a_size && i<nung_tot){

		tol_e(i,a[i],TMath::Cos(b[i]),a,b,a_size,tol);	
		//cout<<"a["<<i<<"] "<<a[i]<<" b["<<i<<"] "<<b[i]<<" tol["<<i<<"] "<<tol[i]<<endl;
		i++;	
	}	
	
	//cout<<" sizof tol "<<sizeof(tol)/sizeof(tol[0]);	
	int tol_size=sizeof(tol)/sizeof(tol[0]);
   	float bin =  tol_least(tol,tol_size);
	//cout<<"bin "<<bin<<endl;
	
	float e_least= tol_least(a,a_size);
	//cout<<"least e "<<e_least<<endl;
	
	float ctheta[nung_tot];
	for(int i=0; i<nung_tot; i++){
		
		ctheta[i]=0;
		//cout<<"ctheta "<<ctheta[i]<<endl;

	
	}
	
	for(int i=0; i<nung_tot; i++){
		
		ctheta[i]=TMath::Cos(b[i]);
		//cout<<"ctheta "<<ctheta[i]<<endl;

	
	}
	TGraph* hctheta = new TGraph(nung_tot,a,ctheta);
	TCanvas* test = new TCanvas("test","",900,900 );

	test->cd();
	//hctheta->SetMarkerStyle(19);
	hctheta->Draw("ap");

	int ctheta_size=sizeof(ctheta)/sizeof(ctheta[0]);
	float ctheta_least= tol_least(ctheta,ctheta_size);
	//cout<<"least ctheta "<<ctheta_least<<endl;
	
  	//int bin_e= (int) 250.0/bin; 
	//int bin_ccthta=(int) 2.0/bin;
        
	//TH2D* hmcp_e_cos_dsig = new TH2D("hmcp_e_cos_dsig","WIMP;cos#theta_{mcp};E_{mcp}",bin_ccthta,-1,1,bin_e,0,250); 
	

	for(int i=0;i<10;i++){

		float dsig_wmp=dsig(a[i],b[i],0);	
		//double bin_cont = hmcp_e_cos_dsig->GetBinContent(bin_no);
		//if(bin_cont==0){
		//	hmcp_e_cos_dsig->Fill(a[i],b[i],dsig_tol);
		//}
	}



	//cout<<"n "<<n<<endl;
	//start including tolerance from here


	
	int nbin=2*nung_tot;
	float edges_a[nbin+1];
	float edges_ctheta[nbin+1];

	//float a_sort[nung_tot], ctheta_sort[nung_tot];

	init_ar(edges_a,nung_tot);
	init_ar(edges_ctheta,nung_tot);
#if 0
	
	for(int i=0; i<sizeof(edges_a)/sizeof(edges_a[0]); i++){
	
		cout<<"edges_a["<<i<<"] "<<edges_a[i]<<endl;
		cout<<"edges_ctheta["<<i<<"] "<<edges_ctheta[i]<<endl;
	
	}
#endif
	//sort(a, a+nung_tot);
	//sort(ctheta, ctheta+nung_tot);

	create_edge(a,tol,edges_a);
	create_edge(ctheta,tol,edges_ctheta);
	sort(edges_a,edges_a+sizeof(edges_a)/sizeof(edges_a[0]));
	sort(edges_ctheta,edges_a+sizeof(edges_a)/sizeof(edges_a[0]));

	for(int i=0; i<sizeof(edges_a)/sizeof(edges_a[0]); i++){
	
		cout<<"edges_a["<<i<<"] "<<edges_a[i]<<endl;
		cout<<"edges_ctheta["<<i<<"] "<<edges_ctheta[i]<<endl;
	
	}
	for(int i=0; i<nung_tot; i++){
	
		cout<<"a["<<i<<"] "<<a[i]<<" tol["<<i<<"] "<<tol[i]<<endl;
		//cout<<"ctheta["<<i<<"] "<<ctheta[i]<<endl;
	
	}
	for(int i=0; i<sizeof(edges_a)/sizeof(edges_a[0]); i++){
	
		//cout<<"edges_a["<<i<<"] "<<edges_a[i]<<endl;
		cout<<"edges_ctheta["<<i<<"] "<<edges_ctheta[i]<<endl;
	
	}
	//int count[nung_tot];
	//for(int i=0; i<nung_tot; i++){
	//
	//	count[i]=0;
	//
	//}

	//TH2F* h = new TH2F("h",";E;cos#theta",nbin,edges_a,nbin,edges_ctheta);
	//TH2F* h1 = new TH2F("h1",";E;cos#theta",nbin,edges_a,nbin,edges_ctheta);
	TH2F* h = new TH2F("h",";E;cos#theta",100,0,250,100,0,1.5);
	TH2F* h1 = new TH2F("h1",";E;cos#theta",100,0,250,10,0,1.5);
	
	//TH1F* h1 = new TH1F("h","",nbin,edges_ctheta);
	//TH1F* h = new TH1F("h","",nbin,edges_ctheta);

    	fill_hist(a,ctheta,tol,h,h1); 
   

	//for(int i=0; i<nung_tot; i++){
	//
	//	h->Fill(a[i]);
	//	h1->Fill(ctheta[i]);
	//
	//}
    
    

 



   

   // for(int ip=0; ip < npfos; ip++){
   //     
   //         if(mcp_pdg[ip]==22 &&  mcp_nparents[ip]==2){
   //         
   //             int m=0;

   //             for(int j= 0; j < mcp_nparents[ip]; j++){
   //             
   //                 if(mcp_pdg[mcp_parentIndex[ip][j]]==11 || mcp_pdg[mcp_parentIndex[ip][j]]==-11) m++; 
   //                 //cout<<"m "<<m<<endl;                
   //             }
   //         
   //             int g=0, v=0, vbar=0;
   //             //if(m==2) {cout<<"m outside loop "<<m<<" in the event ev "<<ev<<endl;}

   //             if(m==2 && mcp_ndaughters[mcp_parentIndex[ip][0]]==3){
   //         
   //                 
   //                 for(int k=0; k< mcp_ndaughters[ip]; k++){
   //                     
   //                      if( mcp_pdg[mcp_daughterIndex[mcp_parentIndex[ip][0]][k]]==22 ) g++;
   //                      if( mcp_pdg[mcp_daughterIndex[mcp_parentIndex[ip][0]][k]]==12 ) v++;
   //                      if( mcp_pdg[mcp_daughterIndex[mcp_parentIndex[ip][0]][k]]==-12 ) vbar++;
   //                             
   //                 }
   //                      
   //                  if(g==1 && v== 1 && vbar==1){ 

   //                      //cout<<"test "<<endl;

   //                      //hmcp_theta->Fill(mcp_theta[ip]);


                            //                      
   //                     // cout<<"e "<<e<<endl;
   //                     // cout<<"theta "<<theta<<endl;
   //                     // cout<<"phi "<<phi<<endl;

                         
   //                 }
   //             }
   //        }        
   //     }


 

	float wt[nung_tot];
	init_ar(wt,nung_tot);
	weight(a,ctheta,tol,wt,nung_tot,b);
	
	//TH2F* hwt = new TH2F("hwt",";E;wt",nbin,edges_a,1000,0,10);
	//TH2F* h2 = new TH2F("h2",";E;wt",nbin,edges_a,1000,0,10);
	TH2F* hwt = new TH2F("hwt",";E;wt",100,0,250,10000,0,0.001);
	TH2F* h2 = new TH2F("h2",";E;wt",100,0,250,10000,0,0.001);

	fill_hist(a,wt,tol,h2,hwt);

	for(int i=0; i<nung_tot; i++){
	
		he_nung->Fill(a[i]);
	
	}
	
	float get_wmp[nung_tot];
	init_ar(get_wmp,nung_tot);
	for(int i=0; i<nung_tot; i++){
	
		int find=he_nung->FindBin(a[i]);
		float get= he_nung->GetBinContent(find);
	
		get_wmp[i]=get*wt[i];

		

		
	
	}

	TGraph* he_wmp = new TGraph(nung_tot,a,get_wmp);
 	TGraphErrors* e_ctheta_wmp = new TGraphErrors(nung_tot,a,ctheta,tol,tol);
	TH1F* he_wmp1 = new TH1F("he_wmp1","",1000,0,250);
	TCanvas* cwmp = new TCanvas("cwmp","E vs cos#theta with tolerance",200,10,700,500);
	cwmp->cd();
	e_ctheta_wmp->SetTitle("E vs Cos#theta including tolerance;E (GeV);Cos#theta");
	e_ctheta_wmp->SetMarkerColor(4);
  	e_ctheta_wmp->SetMarkerStyle(21);
   	e_ctheta_wmp->Draw("ap");

	for(int i=0; i<nung_tot; i++){
	
		he_wmp1->Fill(get_wmp[i]);
	
	}

	//scale he_wmp1 to the pad coordinates
   	Float_t rightmax = 1.1*he_nung->GetMaximum();
   	Float_t scale = gPad->GetUymax()/rightmax;
   	he_nung->SetLineColor(kRed);
   	he_nung->Scale(scale);
	cout<<"Scale "<<scale<<endl;
   	//he_wmp1->Draw("same");

   	
        TCanvas* ce_wmp = new TCanvas("ce_wmp","",900,600);
	
	ce_wmp->cd();
	he_wmp1->SetLineColor(kRed);
	//he_wmp->SetLabelColor(kRed);
	he_nung->Draw("");
	he_wmp1->Draw("same");
	//draw an axis on the right side
   	TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
   	axis->SetLineColor(kRed);
   	axis->SetTextColor(kRed);
   	axis->Draw();






    
        TH1D* hefficiency = (TH1D*) hmcp_pfo_theta->Clone("hefficiency");
        hefficiency->Divide(hmcp_theta);

        TCanvas* cpfo_theta_phi = new TCanvas("phi theta","",600,400);
        TCanvas* cmcp_pfo_theta_phi = new TCanvas("pfo theta","",900,600);
        TCanvas* cefficiency = new TCanvas("eff theta","",900,600);
        TCanvas* cresol = new TCanvas("resol","",900,600);

        TCanvas* cdsig = new TCanvas("dsig","",900,900);
        TCanvas* cwt = new TCanvas("wt","",900,900);
    	//cpfo_theta_phi->Divide(2,2);
    	cresol->Divide(2,2);
    	//cdsig->Divide(1,2);

#if 0
    	cpfo_theta_phi->cd();
    	gPad->SetLogy();
    	gStyle->SetOptFit(1);
    	hmcp_theta->Draw();
        cpfo_theta_phi->SaveAs("cmcp_theta.pdf");

        cmcp_pfo_theta_phi->cd();
    	gPad->SetLogy();
    	gStyle->SetOptFit(1);
    	hmcp_pfo_theta->Draw();
        cmcp_pfo_theta_phi->SaveAs("cmcp_theta_pfo_cond.pdf");

        cefficiency->cd();
    	//gPad->SetLogy();
    	gStyle->SetOptFit(1);
    	hefficiency->Draw();
        cefficiency->SaveAs("efficiency.pdf");

      	cresol->cd(2);
    	//gPad->SetLogy();
    	gStyle->SetOptFit(1);
    	hresol_e->Draw();
       // cpfo_theta_phi->SaveAs("cmcp_theta.pdf");

        cresol->cd(3);
    	//gPad->SetLogy();
    	gStyle->SetOptFit(1);
    	hresol_theta->Draw();

        cresol->cd(4);
    	//gPad->SetLogy();
    	gStyle->SetOptFit(1);
    	hresol_phi->Draw();
        
	cresol->cd(1);
    	//gPad->SetLogy();
    	gStyle->SetOptFit(1);
    	hmcp_pfo_e->Draw();
        cresol->SaveAs("cresol.pdf");
	
	cdsig->cd(1);
	gPad->SetLogy();
	hdsig->Draw();

	cdsig->cd(2);
	hmcp_e_cos->Draw();

#endif
	//cdsig->cd(1);
	//h->Draw();
	//hmcp_e_cos->Draw("COLZ");


#if 0
	cdsig->cd();
	gPad->SetFillColor(33);
	//hmcp_e_cos_dsig->Draw("COLZ");
	h1->Draw("COLZ");
	TH1D * projh2X = h1->ProjectionX();
   	TH1D * projh2Y = h1->ProjectionY();

	TCanvas *c1 = new TCanvas("c1", "c1",900,900);
	c1->Divide(2,1);
        gStyle->SetOptStat(0);

	c1->cd(1);
	projh2X->Draw("COLZ");

	c1->cd(2);
	projh2Y->Draw("COLZ");
#endif
	cdsig->cd();
        gPad->SetFillColor(33);
	h1->SetMarkerStyle(19);
        //hmcp_e_cos_dsig->Draw("COLZ");
        h1->Draw("COLZ");

	cwt->cd();
	cwt->SetLogy();
        gPad->SetFillColor(33);
        //hmcp_e_cos_dsig->Draw("COLZ");
	hwt->SetMarkerStyle(19);
        hwt->Draw("COLZ");



	auto t2 = high_resolution_clock::now();
	auto diff = duration_cast<duration<double>>(t2 - t1); 
	// now elapsed time, in seconds, as a double can be found in diff.count()
	long ms = (long)(1000*diff.count());
	cout<<"time required in ms "<<ms<<endl;
}

float dsig(float e, float cthta, float tol1){
	
	int j0=0,s=500*500,spn=0,mx=120;
	float alpha=1.0/137.0, ke=0.6, sig_an=7.0;//unit of cross sect in picobarn
 
	float dsig_const = (alpha * ke * sig_an * pow(4,j0)*pow((2*spn+1),2))/(16.0 * 3.1416);	
	float x = 2*(e+tol1)/sqrt(s);

	cthta+=tol1;
	float cs=dsig_const*( 1 + pow( ( 1- x ),2) )*pow(1-(4*mx*mx)/((1-x)*s),(0.5+j0))/(x * (1-pow(cthta,2)) );
	
	if(std::isnan(cs)) return 0;
	else return cs;

}


int nung_count(float ener, float cangle, float *a, float *b,int size, float tol1){

	int i=0, count = 0;
	bool m=1, n=1;
	while( i<size && i<nung_tot ){
		n = (TMath::Cos(b[i])<=(cangle+tol1) && TMath::Cos(b[i]) >(cangle-tol1) );
		m = (a[i]<=(ener+tol1) && a[i]>(ener-tol1) );
		if(m && n) count++;
		i++;
		}
	
	return count; 	



}


void tol_e(int nung_num,float e, float cthta,float *a, float *b,int size,  float *tol){

 	float tol1=0.1;	
	float dcs= dsig(e, cthta,0);
	//float dcsN= dsig(e, cthta,0);
	float dcs1= dsig(e, cthta,tol1);
	float dcs2= dsig(e, cthta,-tol1);
//	cout<<"nung_count("<<e<<","<<cthta<<",a, b,"<<tol1<<") "<<nung_count(e, cthta,a, b,tol1)<<endl;
	//cout<<"ctheta "<<cthta<<endl;
	while( (TMath::Abs(dcs-dcs1)) > 0.01 && nung_count(e, cthta,a, b,size,tol1)>=1 &&  (TMath::Abs(dcs-dcs2)) > 0.01 ){
		if(nung_count(e,cthta,a,b,size,tol1)>1) tol1=tol1/10;	
		tol1 = tol1/10;
		dcs1= dsig(e,cthta, tol1);
		dcs2= dsig(e, cthta,-tol1);
		//cout<<"tol1 "<<tol1<<endl;
		//cout<<"dcs-dcs1 "<<(dcs-dcs1)<<" dcs-ds2 "<<(dcs-dcs2)<<" nung_count(e, cthta,a, b,tol1) "<<nung_count(e, cthta,a, b,size,tol1)<<endl;	

	
	}

	tol[nung_num] = tol1;
	//cout<<"nung_num "<<nung_num<<endl;
	//cout<<"tol1 "<<tol1<<endl;
	



}


float tol_least(float *tol, int tol_size){

	//int i=0;
	float temp=1000.0;
	//for(int *it = std::begin(tol); it!=std::end(tol); ++it){
	for(int i = 0; i<tol_size; ++i){
	
		if(temp >= tol[i]){
		       	temp = tol[i];
			//cout<<"i inside tol_least "<<i<<endl;
		}	
		
	
	}
	return temp;

}



void create_edge( float *arr, float *tol, float *edges_a){

#if 0
	for(int ind=0; ind<nung_tot; ind++){	
		edges_a[2*ind]= arr[ind]-tol[ind];
		edges_a[2*ind+1]= arr[ind]+tol[ind];
	}
#endif
	for(int ind=0; ind<nung_tot; ind++){	
		edges_a[2*ind]= arr[ind]-tol[ind];
		edges_a[2*ind+1]= arr[ind]+tol[ind];
	}


}

//void resize(float (&a)[nung_count], float bin){
//
//	while()
//
//
//
//}


void init_ar(float *a, int size){

	for(int i=0; i<size; i++){
	
		a[i]=0.0;
	
	}

}

//void fill_hist(float *a, float *ctheta, float *tol, TH2F* h, TH2F* h1){
void fill_hist(float *a, float *ctheta, float *tol, TH2F* h, TH2F* h1){


	float get;	
	//int size=sizeof(a)/sizeof(a[0]);
	for(int i=0; i<nung_tot; i++){
		
		//if(nung_count())
		//if( a[i]<edges_a[2*i+1] && a[i]>edges_a[2*i]) count[i]++;
		//if(count[i]==0 && a[i]<edges_a[2*i+1] && a[i]>edges_a[2*i]) {
		h->Fill(a[i],ctheta[i]);
	}


	for(int i=0; i<nung_tot; i++){
		//h->Fill(ctheta[i]);
		int find = h->FindBin(a[i],ctheta[i]);
		//int find = h->FindBin(ctheta[i]);
		get= h->GetBinContent(find);

		//cout<<"get before "<<get<<endl;

		if(!(get>2.0)  ){
			//cout<<"dsig("<<a[i]<<","<<ctheta[i]<<")="<<dsig(a[i],ctheta[i],0)*tol[i]*tol[i]*pow(10,5)<<endl;
			//h1->Fill(ctheta[i],a[i], dsig(a[i],ctheta[i],0)*tol[i]*tol[i]*pow(10,5));
			cout<<ctheta[i]<<" "<<a[i]<<" dsig "<<dsig(a[i],ctheta[i],0)<<endl;
			//cout<<ctheta[i]<<endl;
			//h1->Fill(a[i],ctheta[i]);
			h1->Fill(a[i],ctheta[i],dsig(a[i],ctheta[i],0));
			//cout<<"get==1 "<<get<<endl;
		
		
		}
	}

		//	count[i]++;
		//}
	
	
}

void weight(float *a, float *ctheta, float *tol, float *wt,int size, float *b){
	
	float sig_L=28093.1;
	float sig_R =1937.6;
	float sig_tot= 0.5*0.5*sig_L+0.5*0.5*sig_R;
	for(int i=0; i<size; i++){
	
		wt[i]=dsig(a[i],ctheta[i],tol[i])*tol[i]*tol[i]/(nung_count(a[i],ctheta[i],a,b,size,tol[i])*sig_tot/nung_tot);
		cout<<"nung count "<<nung_count(a[i],ctheta[i],a,b,size,tol[i])<<"a["<<i<<"] "<<a[i]<<" ctheta["<<i<<"] "<<ctheta[i]<<endl;
	
	}



}
