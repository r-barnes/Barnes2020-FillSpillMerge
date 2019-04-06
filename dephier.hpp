#ifndef _dephier_hpp_
#define _dephier_hpp_

#include "radix_heap.hpp"
#include <richdem/common/Array2D.hpp>
#include <richdem/common/timer.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <richdem/common/grid_cell.hpp>
#include <richdem/common/constants.hpp>
#include "DisjointDenseIntSet.hpp"
#include "../common/netcdf.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
// #include <unordered_map>
#include <parallel_hashmap/phmap.h>
#include <parallel_hashmap/phmap_utils.h>

namespace richdem::dephier {

//We use a 32-bit integer for labeling depressions. This allows for a maximum of
//2,147,483,647 depressions. This should be enough for most practical purposes.
typedef int32_t dh_label_t;

//Some special valuess
const dh_label_t NO_PARENT = -1;
const dh_label_t NO_VALUE  = -1;

//This class holds information about a depression. Its pit cell and outlet cell
//(in flat-index form) as well as the elevations of these cells. It also notes                                                   //what is flat-index form?
//the depression's parent. The parent of the depression is the outlet through
//which it must flow in order to reach the ocean. If a depression has more than
//one outlet at the same level one of them is arbitrarily chosen; hopefully this                                                  //so, everything should have a parent except for the ocean, right?
//happens only rarely in natural environments.
template<class elev_t>
class Depression {
 public:
  //Flat index of the pit cell, the lowest cell in the depression. If more than
  //one cell shares this lowest elevation, then one is arbitrarily chosen.
  dh_label_t pit_cell = NO_VALUE;
  //Flat index of the outlet cell. If there is more than one outlet cell at this
  //cell's elevation, then one is arbitrarily chosen.
  dh_label_t out_cell = NO_VALUE;
  //Parent depression. If both this depression and its neighbour fill up, this
  //parent depression is the one which will contain the overflow.
  dh_label_t parent   = NO_PARENT;
  //Outlet depression. The metadepression into which this one overflows. Usually
  //its neighbour depression, but sometimes the ocean.
  dh_label_t odep     = NO_VALUE;
  //When a metadepression overflows it does so into the metadepression indicated
  //by `odep`. However, odep must flood from the bottom up. Therefore, we keep
  //track of the `geolink`, which indicates what leaf depression the overflow is
  //initially routed into.
  dh_label_t geolink  = NO_VALUE;
  //Elevation of the pit cell. Since the pit cell has the lowest elevation of
  //any cell in the depression, we initialize this to infinity.
  elev_t  pit_elev = std::numeric_limits<elev_t>::infinity();
  //Elevation of the outlet cell. Since the outlet cell has the lowest elevation
  //of any path leading from a depression, we initialize this to infinity.
  elev_t  out_elev = std::numeric_limits<elev_t>::infinity();
  //The depressions form a binary tree. Each depression has two child
  //depressions: one left and one right.
  dh_label_t lchild = NO_VALUE;
  dh_label_t rchild = NO_VALUE;
  //Indicates whether the parent link is to either the ocean or a depression
  //that links to the ocean
  bool ocean_parent = false;
  //Indicates depressions which link to the ocean through this depression, but
  //are not subdepressions. That is, these ocean-linked depressions may be at
  //the top of high cliffs and spilling into this depression.
  std::vector<dh_label_t> ocean_linked;
  //the label of the depression, for calling it up again
  dh_label_t dep_label = 0;
  //Number of cells contained within the depression and its children
  uint32_t cell_count = 0;
  //Total of elevations within the depression, used in the WLE. Because I think I need to start adding up total elevations before I know the outlet of the depression. 
  //double dep_sum_elevations = 0;
  //Volume of the depression and its children. Used in the Water Level Equation (see below).
  double   dep_vol    = 0;
  //Water currently contained within the depression. Used in the Water Level
  //Equation (see below).
  double   water_vol  = 0;

  //Total elevation of cells contained with the depression and its children
  double total_elevation = 0;
};



//A key part of the algorithm is keeping track of the outlets which connect
//depressions. While each depression has only one outlet, a depression may have
//many inlets. All of the outlets and inlets taken together form a graph which
//we can traverse to determine which way water flows. This class keeps track of
//which cell links two depressions, as well as the elevation of that cell.

//The OutletLink is used as a key for a hashtable which stores information about
//the outlets.
class OutletLink {
 public:
  dh_label_t depa;
  dh_label_t depb;
  OutletLink() = default;
  OutletLink(dh_label_t depa0, dh_label_t depb0) : depa(depa0), depb (depb0) {}
  //This is used to compare two outlets. The outlets are the same regardless of
  //the order in which they store depressions
  bool operator==(const OutletLink &o) const {
    return depa==o.depa && depb==o.depb;
  }
};

//The outlet class stores, again, the depressions being linked as well as
//information about the link
template<class elev_t>
class Outlet {
 public:
  dh_label_t depa;                //Depression A
  dh_label_t depb;                //Depression B
  dh_label_t out_cell = NO_VALUE; //Flat-index of cell at which A and B meet.
  //Elevation of the cell linking A and B
  elev_t  out_elev = std::numeric_limits<elev_t>::infinity();

  Outlet() = default;

  //Standard issue constructor
  Outlet(dh_label_t depa0, dh_label_t depb0, dh_label_t out_cell0, elev_t out_elev0){
    depa       = depa0;
    depb       = depb0;
    if(depa>depb)           //Create a preferred ordering so that comparisons and hashing are faster
      std::swap(depa,depb);
    out_cell   = out_cell0;
    out_elev   = out_elev0;
  }

  //Determines whether one outlet is the same as another. Note that we do not
  //check elevation for this! This is because we'll be using this operator to
  //determine if we've already found an outlet for a depression. We'll look at
  //outlets from lowest to highest, so if an outlet already exists for a
  //depression, it is that depression's lowest outlet.
  bool operator==(const Outlet &o) const {                                                                                              //so beyond just checking, is this somehow preventing it from being recorded if one already exists? How does this work?
    //Outlets are the same if they link two depressions, regardless of the
    //depressions' labels storage order within this class.
    return depa==o.depa && depb==o.depb;
  }
};

//We'll initially keep track of outlets using a hash table. Hash tables require
//that every item they contain be reducible to a numerical "key". We provide
//such a key for outlets here.
template<class elev_t>
struct OutletHash {
  std::size_t operator()(const OutletLink &out) const {
    //Since depa and depb are sorted on construction, we don't have to worry
    //about which order the invoking code called them in and our hash function
    //doesn't need to be symmetric with respect to depa and depb.
    return phmap::HashState().combine(0, out.depa, out.depb);
    //Hash function from: https://stackoverflow.com/a/27952689/752843
    // return out.depa^(out.depb + 0x9e3779b9 + (out.depa << 6) + (out.depa >> 2));    
  }
};



//The regular mod function allows negative numbers to stay negative. This mod
//function wraps negative numbers around. For instance, if a=-1 and n=100, then
//the result is 99.
int ModFloor(int a, int n) {
  return ((a % n) + n) % n;
}


template<class elev_t>
using PriorityQueue = radix_heap::pair_radix_heap<elev_t,uint64_t>;



//Cell is not part of a depression
const dh_label_t NO_DEP = -1; 
//Cell is part of the ocean and a place from which we begin searching for
//depressions.
const dh_label_t OCEAN  = 0;

template<typename elev_t>
using DepressionHierarchy = std::vector<Depression<elev_t>>;

template<class elev_t>
void CalculateMarginalVolumes(DepressionHierarchy<elev_t> &deps, const Array2D<elev_t> &dem, const Array2D<int> &label);

template<class elev_t>
void CalculateTotalVolumes(DepressionHierarchy<elev_t> &deps);



//Calculate the hierarchy of depressions. Takes as input a digital elevation
//model and a set of labels. The labels should have `OCEAN` for cells
//representing the "ocean" (the place to which depressions drain) and `NO_DEP`
//for all other cells.
//
//@param  dem   - 2D array of elevations. May be in any data format.
//
//@return label - A label indiciate which depression the cell belongs to.
//                The indicated label is always the leaf of the depression
//                hierarchy, or the OCEAN.
//
//        flowdirs - A value [0,7] indicated which direction water from the cell
//                   flows in order to go "downhill". All cells have a flow
//                   direction (even flats) except for pit cells.
template<class elev_t, Topology topo>                                                     
DepressionHierarchy<elev_t> GetDepressionHierarchy(
  const Array2D<elev_t> &dem,
  Array2D<dh_label_t>   &label,
  Array2D<int8_t>       &flowdirs
){
  ProgressBar progress;
  Timer timer_overall;
  Timer timer_dephier;
  timer_overall.start();

  timer_dephier.start();

  std::cerr<<"\033[91m##########Getting depression hierarchy\033[39m"<<std::endl;

  //A D4 or D8 topology can be used.
  const int    *dx;
  const int    *dy;
  const int    *dinverse;
  int     neighbours;
  if(topo==Topology::D4){
    dx         = d4x;
    dy         = d4y;
    dinverse   = d4_inverse;
    neighbours = 4;
  } else if(topo==Topology::D8){
    dx         = d8x;
    dy         = d8y;
    dinverse   = d8_inverse;
    neighbours = 8;    
  } else {
    throw std::runtime_error("Unrecognised topology!");
  }

  //Depressions are identified by a number [0,*). The ocean is always
  //"depression" 0. This vector holds the depressions.
  DepressionHierarchy<elev_t> depressions;

  //This keeps track of the outlets we find. Each pair of depressions can only
  //be linked once and the lowest link found between them is the one which is
  //retained.
  typedef phmap::flat_hash_map<OutletLink, Outlet<elev_t>, OutletHash<elev_t>> outletdb_t;
  // typedef std::unordered_map<OutletLink, Outlet<elev_t>, OutletHash<elev_t>> outletdb_t;
  outletdb_t outlet_database;

  //The priority queue ensures that cells are visited in order from lowest to
  //highest. If two or more cells are of equal elevation then the one added last
  //(most recently) is returned from the queue first. This ensures that a single
  //depression gets all the cells within a flat area.
  PriorityQueue<elev_t> pq;

  std::cerr<<"p Adding ocean cells to priority-queue..."<<std::endl;
  //We assume the user has already specified a few ocean cells from which to
  //begin looking for depressions. We add all of these ocean cells to the
  //priority queue now.
  int ocean_cells = 0;
  #pragma omp parallel for collapse(2) reduction(+:ocean_cells)
  for(int y=0;y<dem.height();y++)
  for(int x=0;x<dem.width();x++){
    if(label(x,y)!=OCEAN)
      continue;
    bool has_non_ocean = false;
    for(int n=1;n<=neighbours;n++){
      if(label.inGrid(x+dx[n],y+dy[n]) && label(x+dx[n],y+dy[n])!=OCEAN){
        has_non_ocean = true;
        break;
      }
    }
    if(has_non_ocean){       //If they are ocean cells, put them in the priority queue
      #pragma omp critical
      pq.emplace(dem(x,y), dem.xyToI(x,y));
      ocean_cells++;
    }
  }

  //But maybe the user didn't specify any cells! We'll assume this was mistake
  //and throw an exception. The user can always catch it if they want to.
  if(ocean_cells==0)
    throw std::runtime_error("No initial ocean cells were found!");

  //The 0th depression is the ocean. We add it to the list of depressions now
  //that we're sure there is an ocean!
  { //Use a little scope to avoid having `oceandep` linger around
    auto &oceandep    = depressions.emplace_back();
    //The ocean is deep
    oceandep.pit_elev = -std::numeric_limits<elev_t>::infinity();
    //It's so deep we can't find its bottom
    oceandep.pit_cell = NO_VALUE;
    oceandep.dep_label = 0;
  }


  std::cerr<<"p Finding pit cells..."<<std::endl;

  //Here we find the pit cells of internally-draining regions. We define these
  //to be cells without any downstream neighbours. Note that this means we will
  //identify all flat cells as being pit cells. For DEMs with many flat cells,
  //this will bloat the priortiy queue slightly. If your DEM includes extensive,
  //predictably located flat regions, you may wish to add these in some special
  //way. Alternatively, you could use Barnes (2014, "An Efficient Assignment of
  //Drainage Direction Over Flat Surfaces") as a way of reducing the number of
  //flat cells. Regardless, the algorithm will deal gracefully with the flats it
  //finds and this shouldn't slow things down too much!
  int pit_cell_count = 0;
  progress.start(dem.size());
  #pragma omp parallel for collapse(2) reduction(+:pit_cell_count)
  for(int y=0;y<dem.height();y++)  //Look at all the cells
  for(int x=0;x<dem.width() ;x++){ //Yes, all of them
    ++progress;
    if(label(x,y)==OCEAN)          //Already in priority queue
      continue;
    const auto my_elev = dem(x,y); //Focal cell's elevation
    bool has_lower     = false;    //Pretend we have no lower neighbours
    for(int n=1;n<=neighbours;n++){ //Check out our neighbours
      //Use offset to get neighbour x coordinate, wrapping as needed
      // const int nx = ModFloor(x+dx[n],dem.width()); 
      const int nx = x+dx[n];
      //Use offset to get neighbour y coordinate
      const int ny = y+dy[n];      
      if(!dem.inGrid(nx,ny))  //Is cell outside grid (too far North/South)?
        continue;             //Yup: skip it.
      if(dem(nx,ny)<my_elev){ //Is this neighbour lower than focal cell?
        has_lower = true;     //Make a note of it
        break;                //Don't need to look at additional neighbours
      }
    }
    if(!has_lower){           //The cell can't drain, so it is a pit cell
      pit_cell_count++;       //Add to pit cell count. Parallel safe because of reduction.
      #pragma omp critical    //Only one thread can safely access pq at a time
      pq.emplace(dem(x,y), dem.xyToI(x,y)); //Add cell to pq
    }
  }
  progress.stop();
  std::cerr<<"t Pit cells found in = "<<progress.time_it_took()<<" s"<<std::endl;



  //The priority queue now contains all of the ocean cells as well as all of the
  //pit cells. We will now proceed through the cells by always pulling the cell
  //of lowest elevation from the priority queue. This ensures that, when we find
  //an outlet between two depressions, it is always the lowest outlet linking
  //them.

  //Once two depressions meet, they do not battle for dominance. Rather, we
  //build an invisible, and solely conceptual, wall between them. Each
  //depression continues to grow, by encompassing cells which have not yet been
  //assigned to a depression, until it is surrounded by other depressions on all
  //sides. The depression then contains its pit cell, its outlet, every cell
  //below its outlet, and possibly many cells above its outlet which ultimately
  //drain into its pit cell. (Though note that if the depression has a flat
  //bottom the pit cell is chosen arbitrarily as one of the flat cells.) Later,
  //when we construct the depression hierarchy, we will separate the cells above
  //a depression's outlet into new meta-depressions.

  //In the following we'll temporarily relax our definition of an outlet to mean
  //"the lowest connection between two depressions" rather than "the lowest
  //connection out of a depression". This means depressions may have outlets at
  //many elevations. Later on we'll fix this and some of those outlets will
  //become inlets or the outlets of meta-depressions.

  //The hash table of outlets will dynamically resize as we add elements to it.
  //However, this slows things down a bit. Therefore, we presize the hash set to
  //be equal to be 3x the number of pit cells plus the ocean cell. 3 is just a
  //guess as to how many neighbouring depressions each depression will have. If
  //we get this value too small we lose a little speed due to rehashing. If we
  //get this value too large then we waste space. If we get this value far too
  //large then we run out of RAM.
  outlet_database.reserve(3*(pit_cell_count+1));

  //Visit cells in order of elevation from lowest to highest. If two or more
  //cells are of the same elevation then we visit the one added last (most
  //recently) first.

  std::cerr<<"p Searching for outlets..."<<std::endl;

  progress.start(dem.size());
  while(!pq.empty()){
    ++progress;

    const auto ci    = pq.top_value();     //Copy cell with lowest elevation from priority queue
    const auto celev = pq.top_key();       //Elevation of focal cell
    pq.pop();                              //Remove the copied cell from the priority queue
    auto clabel = label(ci);               //Nominal label of cell
    int cx,cy;
    dem.iToxy(ci,cx,cy);

    if(clabel==OCEAN){
      //This cell is an ocean cell or a cell that flows into the ocean without
      //encountering any depressions on the way. Upon encountering it we do not
      //need to do anything special.
    } else if(clabel==NO_DEP){
      //Since cells label their neighbours and ocean cells are labeled in the
      //initialization, the only way to get to a cell that is still labeled as
      //not being part of a depression is if that cell were added as a pit cell.
      //For each pit cell we find, we make a new depression and label it
      //accordingly. Not all the pit cells originally added will form new
      //depressions as flat cells will relabel their neighbours and the first                                                 
      //cell found in a flat determines the label for the entirety of that flat.
      clabel            = depressions.size();         //In a 0-based indexing system, size is equal to the id of the next flat
      auto &newdep      = depressions.emplace_back(); //Add the next flat (increases size by 1)
      newdep.pit_cell   = dem.xyToI(cx,cy);           //Make a note of the pit cell's location
      newdep.pit_elev   = celev;                      //Make a note of the pit cell's elevation
      newdep.dep_label  = clabel;                     //I am storing the label in the object so that I can find it later and call up the number of cells and volume (better way of doing this?) -- I have since realised I can use the index in the depressions array. So perhaps the label is no longer needed?
      label(ci)         = clabel;                     //Update cell with new label                                                           
    } else {

      //Cell has already been assigned to a depression. In this case, one of two
      //things is true. (1) This cell is on the frontier of our search, in which
      //case the cell has neighbours which have not yet been seen. (2) This cell
      //was part of a flat which has previously been processed by a wavefront
      //beginning at some other cell. In this case, all of this cell's
      //neighbours will have already been seen and added to the priority queue.   
      //However, it is harmless to check on them again.
    }

    //TODO: Update the appropriate depression's cell_count and dep_vol variables                        I did this in the else if above, and then in the if and the else below. I add a cell to the count whenever it is added to the depression and add its elevation to the total elevations.
    //here.                                                                                             Then I calculate the total volume only when we find an outlet (Good way to test this? Print values of volumes only of those that make the outlet queue? I get some negative values sometimes so I may have done something wrong, but what if it's an 'outlet' at the highest point of the depression?)

    //Consider the cell's neighbours
    for(int n=1;n<=neighbours;n++){
      // const int nx = ModFloor(cx+dx[n],dem.width()); //Get neighbour's x-coordinate using an offset and wrapping
      const int nx = cx + dx[n];                      //Get neighbour's y-coordinate using an offset
      const int ny = cy + dy[n];                      //Get neighbour's y-coordinate using an offset
      if(!dem.inGrid(nx,ny))                          //Is this cell in the grid?
        continue;                                     //Nope: out of bounds.
      const auto ni     = dem.xyToI(nx,ny);           //Flat index of neighbour
      const auto nlabel = label(ni);                  //Label of neighbour

      if(nlabel==NO_DEP){                             //Neighbour has not been visited yet 
        label(ni) = clabel;                           //Give the neighbour my label
        pq.emplace(dem(ni), dem.xyToI(nx,ny));        //Add the neighbour to the priority queue
        flowdirs(nx,ny) = dinverse[n];                //Neighbour flows in the direction of this cell
      } else if (nlabel==clabel) {
        //Skip because we are not interested in ourself. That would be vain.
        //Note that this case will come up frequently as we traverse flats since
        //the first cell to be visited in a flat labels all the cells in the
        //flat like itself. All of the other cells in the flat will come off of
        //the priority queue later and, in looking at their neighbours, reach
        //this point.
      } else {
        //We've found a neighbouring depression!

        //Determine whether the focal cell or this neighbour is the outlet of
        //the depression. The outlet is the higher of the two.
        auto out_cell = ci;    //Pretend focal cell is the outlet
        auto out_elev = celev; //Note its height

        if(dem(ni)>out_elev){  //Check to see if we were wrong and the neighbour cell is higher.
          out_cell = ni;       //Neighbour cell was higher. Note it.
          out_elev = dem(ni);  //Note neighbour's elevation
        }

        //We've found an outlet between two depressions. Now we need to
        //determine if it is the lowest outlet between the two.

        //Even though we pull cells off of the priority queue in order of
        //increasing elevation, we can still add a link between depressions that
        //is not as low as it could be. This can happen at saddle points, for
        //instance, consider the cells A-H and their corresponding elevations.
        //Cells in parantheses are in a neighbouring depression
        //     (B) (C) (D)   (256) (197) (329)
        //      A   X   E     228    X    319
        //      H   G   F     255   184   254

        //In this case we are at Cell X. Cells B, C, and D have previously been
        //added. Cell G has added X and X has just been popped. X considers its
        //neighbours in order from A to H. It finds B, at elevation 256, and
        //makes a note that its depression links with B's depression. It then
        //sees Cell C, which is in the same depression as B, and has to update
        //the outlet information between the two depressions.

        const OutletLink olink(clabel,nlabel);      //Create outlet link (order of clabel and nlabel doesn't matter)
        if(outlet_database.count(olink)!=0){        //Determine if the outlet is already present
          auto &outlet = outlet_database.at(olink); //It was. Use the outlet link to get the outlet information
          if(outlet.out_elev>out_elev){             //Is the previously stored link higher than the new one?
            outlet.out_cell = out_cell;             //Yes. So update the link with new outlet cell
            outlet.out_elev = out_elev;             //Also, update the outlet's elevation
          }
        } else {                                    //No preexisting link found; create a new one
          outlet_database[olink] = Outlet<elev_t>(clabel,nlabel,out_cell,out_elev);   
        }
      }

    }
  }
  progress.stop();
  std::cerr<<"t Outlets found in = "<<progress.time_it_took()<<" s"<<std::endl;

  //At this point every cell is associated with the label of a depression. Each
  //depression contains the cells lower than its outlet elevation as well as all
  //cells whose flow ultimately terminates somewhere within the depression. The
  //next order of business is to determine which depressions flow into which
  //other depressions. That is, we need to build a hierarchy of depressions.

  //Since outlets link two depressions, any time we have an outlet we can form a
  //meta-depression whose two children most both fill before the meta-depression
  //itself can spill. This meta-depression has an outlet which differs from
  //either of its children.

  //We can identify the cells belonging to a meta-depression as those which are
  //less than or equal to its outlet elevation, but greater than the elevation
  //of the outlet linking its two children.

  //Since each depression has one and only one parent and at most two children,
  //the depressions will form a binary tree.

  //In order to build the depression hierarchy, it is convenient to visit
  //outlets from lowest to highest.

  //Since the `unordered_set` is, well, unordered. We create a vector into which
  //we will move all of the outlets so we can sort them by elevation. Note that
  //this temporarily doubles the memory required by the program. TODO: Is there
  //a way to avoid this doubling?
  std::vector<Outlet<elev_t>> outlets;

  //Pre-size the vector to avoid expensive copy operations as we expand it
  outlets.reserve(outlet_database.size());

  //Copy the database into the vector
  for(const auto &o: outlet_database)
    outlets.push_back(o.second);

  //It's a little difficult to free memory, but this should do it by replacing
  //the outlet database with an empty database.
  outlet_database = outletdb_t();

  //Sort outlets in order from lowest to highest. Takes O(N log N) time.
  std::sort(outlets.begin(), outlets.end(), [](const Outlet<elev_t> &a, const Outlet<elev_t> &b){
    return a.out_elev<b.out_elev;
  });

  //TODO: For debugging
  for(unsigned int i=0;i<outlets.size()-1;i++)
    assert(outlets.at(i).out_elev<=outlets.at(i+1).out_elev);

  //Now that we have the outlets in order, we'll visit them from lowest to
  //highest. If two outlets are at the same elevation we visit them in an
  //arbitrary order. Each outlet we find is the unique lowest connection between
  //two depressions. We join these depressions to make a meta-depression. The
  //problem is, once we've formed a meta-depression, there may still be many
  //outlets which believe they link to one of the child depressions.

  //To deal with this, we use a Disjoint-Set/Union-Find data structure. This
  //data structure, when passed a depression label as a query, returns the label
  //of the upper-most meta-depression in the chain of parent depressions
  //starting at the query label. The Disjoint-Set data structure has some nice
  //caching properties which, *roughly speaking*, ensure that all queries
  //execute in O(1) time.

  //Presize the DisjointDenseIntSet to twice the number of depressions. Since we
  //are building a binary tree the number of leaf nodes is about equal to the
  //number of non-leaf nodes. The data structure will expand dynamically as
  //needed.
  DisjointDenseIntSet djset(depressions.size());

  std::cerr<<"p Constructing hierarchy from outlets..."<<std::endl;

  //Visit outlets in order of elevation from lowest to highest. If two outlets
  //are at the same elevation, choose one arbitrarily.
  progress.start(outlets.size());
  for(auto &outlet: outlets){
    ++progress;
    auto depa_set = djset.findSet(outlet.depa); //Find the ultimate parent of Depression A
    auto depb_set = djset.findSet(outlet.depb); //Find the ultimate parent of Depression B
    
    //If the depressions are already part of the same meta-depression, then
    //nothing needs to be done.
    if(depa_set==depb_set)
      continue; //Do nothing, move on to the next highest outlet


    if(depa_set==OCEAN || depb_set==OCEAN){
      //If we're here then both depressions cannot link to the ocean, since we
      //would have used `continue` above. Therefore, one and only one of them
      //links to the ocean. We swap them to ensure that `depb` is the one which
      //links to the ocean.
      if(depa_set==OCEAN){
        std::swap(outlet.depa, outlet.depb);
        std::swap(depa_set, depb_set);
      }

      //We now have four values, the Depression A Label, the Depression B Label,
      //the Depression A MetaLabel, and the Depression B MetaLabel. We know that
      //the Depression B MetaLabel is OCEAN. Depression B Label is the label of
      //the actual depression this outlet links to, not the meta-depressions of
      //which it is a part. Depression A MetaLabel is the meta-depression that
      //has just found a path to the ocean via Depression B. Depression A Label
      //is some value we don't care about.

      //What we will do is link Depression A MetaLabel to Depression B.
      //Depression B ultimately terminates in the ocean, but the only way to get
      //there in real-life is to crawl into Depression B, not into its meta-
      //depression. At this point its meta-depression is the ocean, so crawling
      //into the meta-depression would form a direct link to the ocean, which is
      //not realistic.

      //Get a reference to Depression A MetaLabel.
      auto &dep = depressions.at(depa_set);

      //If this depression has already found the ocean then don't merge it
      //again. (TODO: Richard)
      // if(dep.out_cell==OCEAN)
        // continue;

      //Ensure we don't modify depressions that have already found their paths
      assert(dep.out_cell==-1);
      assert(dep.odep==NO_VALUE);            

      //Point this depression to the ocean through Depression B Label
      dep.parent       = outlet.depb;        //Set Depression Meta(A) parent
      dep.out_elev     = outlet.out_elev;    //Set Depression Meta(A) outlet elevation                                     
      dep.out_cell     = outlet.out_cell;    //Set Depression Meta(A) outlet cell index
      dep.odep         = outlet.depb;        //Depression Meta(A) overflows into Depression B
      dep.ocean_parent = true;
      dep.geolink      = outlet.depb;        //Metadepression(A) overflows, geographically, into Depression B
      depressions.at(outlet.depb).ocean_linked.emplace_back(depa_set);
      djset.mergeAintoB(depa_set,OCEAN); //Make a note that Depression A MetaLabel has a path to the ocean
    } else {
      //Neither depression has found the ocean, so we merge the two depressions
      //into a new depression.
      auto &depa          = depressions.at(depa_set); //Reference to Depression A MetaLabel
      auto &depb          = depressions.at(depb_set); //Reference to Depression B MetaLabel

      //Ensure we haven't already given these depressions outlet information
      assert(depa.odep==NO_VALUE);     
      assert(depb.odep==NO_VALUE);

      const auto newlabel = depressions.size();       //Label of A and B's new parent depression
      depa.parent   = newlabel;        //Set Meta(A)'s parent to be the new meta-depression
      depb.parent   = newlabel;        //Set Meta(B)'s parent to be the new meta-depression
      depa.out_cell = outlet.out_cell; //Note that this is Meta(A)'s outlet
      depb.out_cell = outlet.out_cell; //Note that this is Meta(B)'s outlet
      depa.out_elev = outlet.out_elev; //Note that this is Meta(A)'s outlet's elevation
      depb.out_elev = outlet.out_elev; //Note that this is Meta(B)'s outlet's elevation
      depa.odep     = depb_set;        //Note that Meta(A) overflows, logically, into Meta(B)
      depb.odep     = depa_set;        //Note that Meta(B) overflows, logically, into Meta(A)
      depa.geolink  = outlet.depb;     //Meta(A) overflows, geographically, into B
      depb.geolink  = outlet.depa;     //Meta(B) overflows, geographically, into A
   
      //Be sure that this happens AFTER we are done using the `depa` and `depb`
      //references since they will be invalidated if `depressions` has to
      //resize!
      const auto depa_pitcell_temp = depa.pit_cell;

      auto &newdep     = depressions.emplace_back();                                                                       //is it right to create a new depression for the metadepression like this?
      newdep.lchild    = depa_set;
      newdep.rchild    = depb_set; 
      newdep.dep_label = newlabel;
      newdep.pit_cell  = depa_pitcell_temp;

      djset.mergeAintoB(depa_set, newlabel); //A has a parent now
      djset.mergeAintoB(depb_set, newlabel); //B has a parent now
    }
  }
  progress.stop();

  std::cerr<<"t Time to construct Depression Hierarchy = "<<timer_dephier.stop()<<" s"<<std::endl;


  //At this point we have a 2D array in which each cell is labeled. This label
  //corresponds to either the root node (the ocean) or a leaf node of a binary
  //tree representing the hierarchy of depressions.

  //The labels array has been modified in place. The depression hierarchy is
  //returned.

  Timer timer_volumes;
  timer_volumes.start();

  CalculateMarginalVolumes(depressions, dem, label);

  // { //Depression filling code
  //   Array2D<elev_t> dhfilled(dem);
  //   for(int i=0;i<(int)dem.size();i++)
  //     dhfilled(i) = dem(i);

  //   //Get the marginal depression cell counts and total elevations
  //   progress.start(dem.size());
  //   for(unsigned int i=0;i<dem.size();i++){
  //     ++progress;
  //     auto clabel        = label(i);
      
  //     if(clabel==OCEAN)
  //       continue;

  //     while(depressions[clabel].parent!=OCEAN && !depressions[clabel].ocean_parent)
  //       clabel = depressions[clabel].parent;

  //     if(dem(i)<depressions[clabel].out_elev)
  //       dhfilled(i) = depressions[clabel].out_elev;
  //   }
  //   progress.stop();

  //   SaveAsNetCDF(dhfilled,"/z/out-dhfilled.nc","value");
  // }

  CalculateTotalVolumes(depressions);

  std::cerr<<"t Time to calculate volumes = "<<timer_volumes.stop()<<" s"<<std::endl;
  std::cerr<<"t Total time in depression hierarchy calculations = "<<timer_overall.stop()<<" s"<<std::endl;

  return depressions;
}



template<class elev_t>
void CalculateMarginalVolumes(
  DepressionHierarchy<elev_t> &deps,
  const Array2D<elev_t>   &dem,
  const Array2D<int>      &label
){
  ProgressBar progress;

  std::cerr<<"p Calculating depression marginal volumes..."<<std::endl;

  //Get the marginal depression cell counts and total elevations
  progress.start(dem.size());
  #pragma omp parallel default(none) shared(progress,deps,dem,label)
  {
    std::vector<uint32_t> cell_counts     (deps.size(), 0);
    std::vector<double>   total_elevations(deps.size(), 0);

    #pragma omp for
    for(unsigned int i=0;i<dem.size();i++){
      ++progress;
      const auto my_elev = dem(i);
      auto clabel        = label(i);
      
      while(clabel!=OCEAN && my_elev>deps.at(clabel).out_elev)
        clabel = deps[clabel].parent;

      if(clabel==OCEAN)
        continue;

      cell_counts[clabel]++;
      total_elevations[clabel] += dem(i);
    }

    #pragma omp critical
    for(unsigned int i=0;i<deps.size();i++){
      deps[i].cell_count      += cell_counts[i];
      deps[i].total_elevation += total_elevations[i];
    }
  }
  progress.stop();
}



template<class elev_t>
void CalculateTotalVolumes(
  DepressionHierarchy<elev_t> &deps
){
  ProgressBar progress;

  std::cerr<<"p Calculating depression total volumes..."<<std::endl;
  //Calculate total depression volumes and cell counts
  progress.start(deps.size());
  for(int d=0;d<(int)deps.size();d++){
    ++progress;

    auto &dep = deps.at(d);
    if(dep.lchild!=NO_VALUE){
      assert(dep.rchild!=NO_VALUE); //Either no children or two children
      assert(dep.lchild<d);         //ID of child must be smaller than parent's
      assert(dep.rchild<d);         //ID of child must be smaller than parent's
      dep.cell_count      += deps.at(dep.lchild).cell_count;
      dep.total_elevation += deps.at(dep.lchild).total_elevation;
      dep.cell_count      += deps.at(dep.rchild).cell_count;
      dep.total_elevation += deps.at(dep.rchild).total_elevation;
    }
    //This has to be after the foregoing because the cells added by the if-
    //clauses have additional volume above their spill elevations that cannot be
    //counted simply by adding their volumes to their parent depression.
    dep.dep_vol = dep.cell_count*static_cast<double>(dep.out_elev)-dep.total_elevation;

    assert(dep.lchild==NO_VALUE || deps.at(dep.lchild).dep_vol+deps.at(dep.rchild).dep_vol<=dep.dep_vol);
  }
  progress.stop();
}



//Utility function for doing various relabelings based on the depression
//hierarchy.
template<class elev_t>
void LastLayer(Array2D<dh_label_t> &label, const Array2D<float> &dem, const DepressionHierarchy<elev_t> &depressions){
  #pragma omp parallel for collapse(2)
  for(int y=0;y<label.height();y++)
  for(int x=0;x<label.width();x++){
    auto mylabel = label(x,y);
    while(true){
      if(dem(x,y)>=depressions.at(mylabel).out_elev)
        mylabel = depressions.at(mylabel).parent;
      else {
        if(mylabel!=0)
          mylabel = -3;
        break;
      }
    }
    label(x,y) = mylabel;
  }
}

}

#endif

