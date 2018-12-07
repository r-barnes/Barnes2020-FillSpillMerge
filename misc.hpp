#ifndef _dephier_misc_
#define _dephier_misc_

//Utility function for doing various relabelings based on the depression
//hierarchy.
template<class elev_t>
void FillDepressions(rd::Array2D<label_t> &label, rd::Array2D<float> &dem, const std::vector<Depression<elev_t>> &depressions){
  #pragma omp parallel for collapse(2)
  for(int y=0;y<label.height();y++)
  for(int x=0;x<label.width();x++){
    auto mylabel = label(x,y);
    auto myelev  = dem(x,y);
    while(true){
      if(dem(x,y)>=depressions.at(mylabel).out_elev){
        mylabel = depressions.at(mylabel).parent;
        myelev  = depressions.at(mylabel).out_elev;
      } else {
        break;
      }
    }
    label(x,y) = mylabel;
  }
}

#endif
