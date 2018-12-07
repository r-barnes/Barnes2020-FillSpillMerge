#include "Array2D.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>

constexpr double SQ2 = std::sqrt(2.0);

//1 2 3
//0   4
//7 6 5
//                      0  1  2  3 4 5 6  7
const int dx8[8]       = {-1,-1, 0, 1,1,1,0,-1};
const int dy8[8]       = {0, -1,-1,-1,0,1,1, 1};
const double dr8[8]    = {1,SQ2,1,SQ2,1,SQ2,1,SQ2};
const int d8inverse[8] = {4,  5, 6, 7,0,1,2, 3};

//  1
//0   2
//  3
//                      0  1 2 3
const int dx4[8]       = {-1, 0,1,0};
const int dy4[8]       = { 0,-1,0,1};
const double dr4[4]    = { 1, 1,1,1};
const int d4inverse[4] = { 2, 3,0,1};

const int    *const dx       = dx8;
const int    *const dy       = dy8;
const int    *const dinverse = d8inverse;
const double *const dr       = dr8;
const int neighbours         = 8;


const float  OCEAN_LEVEL = 0;
const int8_t NO_FLOW     = -1;

class GridCell {
 public:
  int x, y;
  GridCell(const int x0, const int y0){
    x = x0;
    y = y0;
  }
};


class GridCellZ {
 public:
  int x, y;
  float z;
  GridCellZ(const int x0, const int y0, const float z0){
    x = x0;
    y = y0;
    z = z0;
  }
  bool operator>(const GridCellZ& a) const {
    return z>a.z; //Less than sorts the queue in reverse
  }
};



void ProcessDepression(
  int cx0,
  int cy0,
  Array2D<bool>        &processed,
  const Array2D<float> &topo,
  const Array2D<int>   &flowdirs,
  Array2D<float> &wtd,
  Array2D<float> &hydro_surf,
  int outlet
){
  std::vector<int> depression_cells;

  if(outlet==-1){
    int cx = cx0;
    int cy = cy0;
    while(true){
      const int n  = flowdirs(cx,cy);
      const int nx = cx+dx[n];
      const int ny = cy+dy[n];
      //If we start going back downhill, we've reached the outlet
      if(topo(nx,ny)<topo(cx,cy)){
        outlet = topo.xyToI(cx,cy);
        break;
      }
    }
  }

  const auto outlet_elev = topo(outlet);

  std::priority_queue<GridCellZ, std::vector<GridCellZ>, std::greater<GridCellZ> > pq;

  const int original_cx0 = cx0;
  const int original_cy0 = cy0;

  //Climb into depression
  pq.emplace(cx0, cy0, hydro_surf(c0));
  while(!pq.empty()){
    const auto c = pq.top();

    bool has_lower = false;
    for(int n=0;n<neighbours;n++){
      const int nx = c.x+dx[n];
      const int ny = c.y+dy[n];
      if(hydro_surf(nx,ny)<hydro_surf(c.x,c.y))
        has_lower = true;
      pq.emplace(nx,ny,hydro_surf(nx,ny));
    }
    if(!has_lower){
      cx0 = c.x;
      cy0 = c.y;
      c0  = topo.xyToI(c.x,c.y);
      break;
    }
  }

  //Climb out of depression along entry path
  {
    int cx = cx0;
    int cy = cy0;
    while(true){
      if(cx==original_cx0 && cy==original_cy0)
        break;
      depression_cells.push_back(topo.xyToI(cx,cy));
      const auto n = flowdirs(cx,cy);
      cx           = cx+dx[n];
      cy           = cy+dy[n];
    }
  }

  while(!depression_cells.empty()){
    if(hydro_surf(c)>topo(c))
      break;
    if(wtd(c)==0) //We ran out of water
      return;
    if(!depression_cells.empty()){
      const int n = depression_cells.back();
      depression_cells.pop_back();
      wtd(n) += wtd(c);
      wtd(c)  = 0;
    }
  }


  double water_volume = wtd(cx0,cy0);
  wtd(cx0,cy0) = 0;



  //Climb out of depression
  pq.emplace(cx0, cy0, hydro_surf(c0));
  double volume         = 0;
  double cell_area      = 0;
  double base_area      = 0;
  double previous_level = 0;
  while(!pq.empty()){
    const auto c = pq.top();

    if(wtd(c.x,c.y)>0){
      water_volume += wtd(c.x,c.y); //Gather water, especially on flats.
      wtd(c.x,c.y)  = 0;
    }

    volume        += cell_area*(hydro_surf(c.x,c.y)-previous_level);
    previous_level = hydro_surf(c.x,c.y);
    cell_area     += 1;
    base_area     += 1*hydro_surf(c.x,c.y); //Cell area * Cell Elevation
    // WaterVolume=(WaterEl-Aelev)*Aarea+(WaterEl-Belev)*Barea+(WaterEl-Celev)*Carea+...
    // 0=(WaterEl-Aelev)*Aarea+(WaterEl-Belev)*Barea+(WaterEl-Celev)*Carea+...-WaterVolume
    // 0=(Aarea+Barea+Carea+...)*WaterEl-(Aarea*Aelev+Barea*Belev+Carea*Celev+...)-WaterVolume
    // 0=TotalArea*WaterEl-BaseAreaTimesBaseHeight-WaterVolume
    // WaterEl=(BaseAreaTimesBaseHeight+WaterVolume)/TotalArea
    if(volume>water_volume){
      const auto water_elev = (base_area+water_volume)/cell_area;
      for(const auto &dc: depression_cells)
        hydro_surf(dc) = water_elev;
      return;
    }

    depression_cells.push_back(topo.xyToI(c.x,c.y));

    for(int n=0;n<neighbours;n++){
      const int nx = c.x+dx[n];
      const int ny = c.y+dy[n];
      if(processed(nx,ny))
        continue;
      if(hydro_surf(nx,ny)<hydro_surf(c.x,c.y)){
        outlet = true;
        break;
      }
      processed(nx,ny) = true;
      pq.emplace(nx,ny,hydro_surf(nx,ny));
    }

    if(outlet){
      const auto outlet_elev = hydro_surf(c.x,c.y);
      for(const auto &dc: depression_cells)
        hydro_surf(dc) = outlet_elev;
      water_volume -= volume;

      if(hydro_surf(c.x,c.y)<outlet_elev)
        ProcessDepression(c.x,c.y,processed,topo,wtd,hydro_surf);
      else 
        wtd(outlet) = water_volume;
      return;
    }
  }
}



void SurfaceWater(const Array2D<float> &topo, Array2D<float> &wtd){
  //Floating-point math means that adding and subtracting water from the `wtd`
  //offset surface will yield non-flat lake surfaces. Therefore, we create a
  //hydraulic surface which captures the level of the water.
  Array2D<float> hydro_surf(topo.width(),topo.height(),0);
  for(int i=0;i<topo.size();i++)
    hydro_surf(i) = topo(i);

  //Our first step is to move all of the water downstream into pit cells. To do
  //so, we calculate steepest-descent flow directions for each cell (though we
  //could also use MFD). Later, we will calculate a new set of flow directions
  //which provide paths out of depressions. Using such flow directions now would
  //mean that some water on the slopes of depressions would not automatically
  //drain to the middle complicating the algorithm later.

  //The flowdirs array points downstream. We use int8_t to save memory
  Array2D<int8_t>  flowdirs(topo.width(),topo.height(),NO_FLOW);

  //Find the steepest descent neighbour for each cell
  #pragma omp parallel for collapse(2)
  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width();x++){
    double greatest_slope = 0;
    auto   max_n          = NO_FLOW;
    for(int n=0;n<neighbours;n++){
      const int nx = x+dx[n];
      const int ny = y+dy[n];
      const auto slope = (topo(x,y)-topo(nx,ny))/dr[n];
      if(slope>greatest_slope){
        greatest_slope = slope;
        max_n          = n;
      }
    }
    flowdirs(x,y) = max_n;
  }

  //Calculate how many upstream cells flow into each cell
  Array2D<char>  dependencies(topo.width(),topo.height(),0);
  #pragma omp parallel for collapse(2)
  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width();x++){
    if(flowdirs(x,y)==NO_FLOW)
      continue;
    for(int n=0;n<neighbours;n++){
      const int nx = x+dx[n];
      const int ny = y+dy[n];
      if(flowdirs(nx,ny)==dinverse[n])
        dependencies(x,y)++;
    }
  }

  //Find the peaks. These are the cells into which no other cells pass flow. We
  //know the flow accumulation of the peaks without having to perform any
  //recursive calculations; they just pass flow downstream. From the peaks, we
  //can begin a breadth-first traversal in the downstream direction by adding
  //each cell to the frontier/queue as its dependency count drops to zero.
  std::queue<int> q;
  #pragma omp parallel for //TODO: Use a custom reduction
  for(int i=0;i<topo.size();i++)
    if(dependencies(i)==0 && flowdirs(i)!=-1)  //Is it a peak?
      #pragma omp critical
      q.emplace(i);         //Yes.

  //Starting with the peaks, pass flow downstream
  while(!q.empty()){
    const auto c = q.front();          //Copy focal cell from queue
    q.pop();                           //Clear focal cell from queue

    //Coordinates of downstream neighbour, if any
    const auto n  = flowdirs(c); 

    //If downstream neighbour is the ocean, we drop our water into it and the
    //ocean is unaffected.
    if(wtd(c)<=OCEAN_LEVEL){
      wtd(c) = 0;
      continue;
    }

    //If we have water, pass it downstream.
    if(wtd(c)>0){
      wtd(n) += wtd(c);
      wtd(c)  = 0;
    }

    //Decrement the neighbour's dependencies. If there are no more dependencies,
    //we can process the neighbour.
    if(--dependencies(n)==0)
      q.emplace(n);                   //Add neighbour to the queue
  }


  //At this point the water is located at the pit cells of all of the
  //depressions. Now we'll use the Priority-Flood method to determine a
  //processing order for depressions.

  //Indicates whether a cell has previously been processed by the algorithm.
  //This prevents the algorithm from going in an infinite loop.
  Array2D<bool> processed(topo.width(),topo.height(),false);

  //The priority queue ensures that cells are visited in order from lowest to
  //highest
  std::priority_queue<GridCellZ, std::vector<GridCellZ>, std::greater<GridCellZ> > pq;

  ////////////////////
  //Locate depressions
  ////////////////////

  //We start by adding all the edge cells to the priority queue
  for(int y=0;y<topo.height();y++){
    pq.emplace(0,             y,topo(0,             y));
    pq.emplace(topo.width()-1,y,topo(topo.width()-1,y));
  }
  for(int x=0;x<topo.width();x++){
    pq.emplace(x, 0,               topo(x,               0));
    pq.emplace(x, topo.height()-1, topo(x, topo.height()-1));
  }

  //Visit cells in order from lowest to highest generating flow directions as we
  //go.
  while(!pq.empty()){
    const auto c = pq.top();       //Copy cell with lowest elevation from priority queue
    pq.pop();                      //Remove the copied cell from the priority queue

    for(int n=0;n<neighbours;n++){
      const int nx = c.x + dx[n];  //Use focal cell's x-coordinate plus offset to get neighbour x-coordinate
      const int ny = c.y + dy[n];  //Use focal cell's y-coordinate plus offset to get neighbour y-coordinate
      if(!topo.inGrid(nx,ny))      //Neighbour cell is out of bounds
        continue;
      if(processed(nx,ny))         //Cell has already been visited
        continue;

      processed(nx,ny) = true;     //Mark cell as having been visited

      //Regardless of whether it was in a depression or not, we make a note that
      //the neighbour cell's "downstream" flow direction was this focal cell.
      //Roughly speaking, for cells which are not in depressions, this
      //corresponds to flow directions similar to D8. Cells in depressions will
      //drain to the depression's deepest point (the pit cell) and the deepest
      //point will drain via the steepest path to the depression's spill point.
      flowdirs(nx,ny) = n;

      pq.emplace(nx,ny,topo(nx,ny));
    }
  }

  //We now have flow directions which ensure that each cell has a path to the
  //edge of the DEM. We now process the cells from the highest points down, thus
  //ensuring that water is passed correctly from upstream to downstream lakes.

  ////////////////////
  //Flow Accumulation
  ////////////////////

  //Calculate how many upstream cells flow into each cell
  dependencies.setAll(0);
  #pragma omp parallel for collapse(2)
  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width();x++){
    if(flowdirs(x,y)==NO_FLOW)
      continue;
    for(int n=0;n<neighbours;n++){
      const int nx = x+dx[n];
      const int ny = y+dy[n];
      if(flowdirs(nx,ny)==dinverse[n])
        dependencies(x,y)++;
    }
  }

  //Find the peaks. These are the cells into which no other cells pass flow. We
  //know the flow accumulation of the peaks without having to perform any
  //recursive calculations; they just pass flow downstream. From the peaks, we
  //can begin a breadth-first traversal in the downstream direction by adding
  //each cell to the frontier/queue as its dependency count drops to zero.

  //Note that `q` is empty, since that was a condition to end a `while` loop
  //above.
  #pragma omp parallel for //TODO: Use a custom reduction
  for(int i=0;i<topo.size();i++)
    if(dependencies(i)==0 && flowdirs(i)!=-1)  //Is it a peak?
      #pragma omp critical
      q.emplace(i);         //Yes.


  //Starting with the peaks, pass flow downstream
  while(!q.empty()){
    const auto c = q.front();         //Copy focal cell from queue
    q.pop();                          //Clear focal cell from queue

    const auto n = flowdirs(c);       //Find the "downstream" neighbour of the focal cell

    if(hydro_surf(n)>hydro_surf(c) && wtd(c)!=0)  //If the neighbour is higher than we are
      ProcessDepression(c,processed,topo,wtd,hydro_surf);

    if(wtd(c)>=0){
      wtd(n) += wtd(c);
      wtd(c)  = 0;
    }

    //Decrement the neighbour's dependencies. If there are no more dependencies,
    //we can process the neighbour.
    if(--dependencies(n)==0)
      q.emplace(n);                   //Add neighbour to the queue
  }
}






int main(int argc, char **argv){
  if(argc!=3){
    std::cout<<"Syntax: "<<argv[0]<<" <Input file directory> <Run Type>"<<std::endl;
    return -1;
  }

  std::string dir      = argv[1];
  std::string run_type = argv[2];

  Array2D<float> rech  (dir+"/Mad_020500_rech_rotated.nc",   "value");   //Recharge (Percipitation minus Evapotranspiration)
  Array2D<float> temp  (dir+"/Mad_020500_temp_rotated.nc",   "value");   //Air temperature   - Used with fslope to generate efolding depth
  Array2D<float> fslope(dir+"/Mad_020500_fslope_rotated.nc", "value");   //100/(1+150*slope) - Used with fslope to generate efolding depth
  Array2D<float> topo  (dir+"/Mad_020500_topo_rotated.nc",   "value");   //Terrain height
  Array2D<float> ksat  (dir+"/Mad_ksat_rotated.nc",          "value");   //Hydrologic conductivity
  Array2D<float> wtd   (dir+"/Mad_021000_wtd_rotated.nc",    "value");

  if(run_type=="equilibrium")
    throw std::runtime_error("equilibrium not implemented!");
  else if (run_type=="transient"){
    //Pass
  } else
    throw std::runtime_error("Expected 'equilibrium' or 'transient'!");

  //TODO: Make sure all files have same dimensions

  for(int i=0;i<wtd.size();i++){
    wtd(i) = 0;
    if(topo(i)>0)
      wtd(i) = 1;
  }

  SurfaceWater(topo, wtd);

  std::ofstream fout("/z/out.dat");
  fout.write(reinterpret_cast<const char*>(wtd.data), wtd.size()*sizeof(float));

  return 0;
}