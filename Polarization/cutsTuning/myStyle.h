
void SetMyStyle(){
	//TCanvas *c1 = new TCanvas("c1","c1",50,50,600,600);
	//c1->SetFillColor(10);
	gStyle->SetOptFit(1);             //1= display fitted result, 0 = no
	gStyle->SetOptStat(1);            //1= display the statistics, 0 = no
	gStyle->SetTitleFont(22);         // Font for X, Y, and Top title
	gStyle->SetStatFont(22);          // Font for statistics
	gStyle->SetLabelFont(22,"X");     // Font for X label
	gStyle->SetLabelFont(22,"Y");     // Font for Y label
	gStyle->SetTitleXOffset(1.0);     // X title offset
	gStyle->SetTitleYOffset(1.0);     // Y title offset
	gStyle->SetHistLineWidth(2);      // Histogram line width
	gStyle->SetStatX(0.9);            // Centroid X position of statistics box
	gStyle->SetStatY(0.9);            // Centroid Y position of statistics box
	gStyle->SetTitleX(0.15);          // Centroid X position of title box
	gStyle->SetTitleY(0.96);          // Centroid Y position of title box

	gPad->SetLeftMargin(0.15);        // pad left margin  for writing Y title
	gPad->SetBottomMargin(0.15);      // pad bottom margin for writing X title
	gPad->SetRightMargin(0.15);       // pad right margin  for writing Y title
	gPad->SetTopMargin(0.15);         // pad top margin for writing X title
}
