
void Fill1DHist(TH1D *hist, const double &value){
  if (value<hist->GetXaxis()->GetXmin()){
    hist->Fill(hist->GetXaxis()->GetXmin());
  }
  else if (value<hist->GetXaxis()->GetXmax()){
    hist->Fill(value);
  }
  else{
    hist->Fill(hist->GetXaxis()->GetXmax()-0.000001);
  }
}


void Fill1DWHist(TH1D *hist, const double &value, const double &wei){
  if (value<hist->GetXaxis()->GetXmin()){
    hist->Fill(hist->GetXaxis()->GetXmin(), wei);
  }
  else if (value<hist->GetXaxis()->GetXmax()){
    hist->Fill(value, wei);
  }
  else{
    hist->Fill(hist->GetXaxis()->GetXmax()-0.000001, wei);
  }
}

