#ifndef _df_flow_hpp_
#define _df_flow_hpp_

#include "dephier.hpp"
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
#include <richdem/common/Array2D.hpp>
#include <richdem/common/timer.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <richdem/depressions/depressions.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace richdem {

namespace dephier {

const int *const dx         = d8x;
const int *const dy         = d8y;
const int *const dinverse   = d8_inverse;
const int        neighbours = 8;

const double FP_ERROR = 1e-4;

const float  OCEAN_LEVEL = 0;  //ocean_level in the topo file must be lower than any non-ocean cell. 

//TODO: Read this
//We need an efficient way to determine the volume of a depression. To do so, we
//note that if the outlet of a depression containing cells of elevations
//$\{a,b,c,d\}$ is known and at elevation $o$, then the volume of the depression
//is $(o-a)+(o-b)+(o-c)+(o-d)=4o-a-b-c-d=4o-sum(elevations)$. This says that, if
//we keep track of the number of cells in a depression and their total
//elevation, it is possible to calculate the volume of a depression at any time
//based on a hypothetical outlet level. We call this the Water Level Equation.                                                          

//Our strategy will be to keep track of the necessary components of the water
//level equation for each depression and each outlet. Then as the outlets of
//depression become known we can use the WLE to calculate their volumes. The
//volume of a depression is its WLE minus the WLEs of its child depressions.


  //TODO: 0. Calculate the number of cells within and volume of each depression.   DONE 
  //This will take place inside of GetDepressionHierarchy.
 

  //TODO 1. DONE. Get flow directions for all cells. 

  //TODO 2. Perform a flow accumulation moving water downhill and filling             NEARLY DONE - think about oceans
  //groundwater as you go. Look at `misc/main.cpp`. Recall that we did this by
  //counting depressions, finding peaks, and using a queue to control a breadth-
  //first traversal from the peaks downwards.

  //TODO 3. As part of the above, when flow can't go downhill any farther, add        DONE
  //it to the `water_vol` for the appropriate depression. Use the labels array
  //to determine the appropriate depression.

  //TODO 4/5. Perform a depth-first post-order traversal of the depression
  //hierarchy (start with depressions for which `parent==NO_PARENT`. When you
  //reach the leaves if `water_vol>dep_vol` then try to overflow into the
  //neighbouring depression if its `water_vol<dep_vol`. To do so search
  //neighbour cells of `out_cell` for the lowest cell labeled `odep` and follow
  //that one's flow path until it terminates. Add excess water to that
  //depression's `water_vol`. After both child depressions have been visited
  //they will be finished trying to share their water and their excess water is
  //added to their parent's `water_vol`. Repeat the overflow attempt.

  //TODO 6: Adjust the hydrologic elevations. If a depression has `water_vol>0`
  //then all of its child depression are full. What remains is to use a
  //priority-queue and the Water Level Equation (see dephier.cpp) to find which
  //cells should be flooded. If a depression is a leaf depression (no children)
  //then start the PQ at the pit_cell. Otherwise, start at the pit cell of any
  //child depression since all children are flooded to at least the level of the
  //`out_elev` connecting the uppermost two children in the hierarchy. Cells
  //should be added to the queue with elevation equal to at least the `out_elev`
  //of the depression's children and, if `water_vol<dep_vol` should not exceed
  //the `out_elev` of the depression itself.



template<class elev_t, class wtd_t>
void MoveWaterIntoPits(
  const rd::Array2D<elev_t>    &topo,
  rd::Array2D<wtd_t>           &wtd,
  const rd::Array2D<int>       &label,
  DepressionHierarchy<elev_t>  &deps,
  const rd::Array2D<flowdir_t> &flowdirs
){
  rd::Timer timer;
  rd::ProgressBar progress;
  timer.start();

  //Our first step is to move all of the water downstream into pit cells. To do
  //so, we use the steepest-descent flow directions provided by the depression
  //hierarchy code

  std::cerr<<"p Moving surface water downstream..."<<std::endl;

  //Calculate how many upstream cells flow into each cell
  rd::Array2D<char>  dependencies(topo.width(),topo.height(),0);
  #pragma omp parallel for collapse(2)
  for(int y=0;y<topo.height();y++)
  for(int x=0;x<topo.width(); x++)
  for(int n=1;n<=neighbours;n++){     //Loop through neighbours
    const int nx = x+dx[n];           //Identify coordinates of neighbour
    const int ny = y+dy[n];
    if(!topo.inGrid(nx,ny))
      continue;    
    if(flowdirs(nx,ny)==dinverse[n])  //Does my neighbour flow into me?
      dependencies(x,y)++;            //Increment my dependencies
  }

  //Find the peaks. These are the cells into which no other cells pass flow (i.e. 0 dependencies). We
  //know the flow accumulation of the peaks without having to perform any
  //recursive calculations; they just pass flow downstream. From the peaks, we
  //can begin a breadth-first traversal in the downstream direction by adding
  //each cell to the frontier/queue as its dependency count drops to zero.
  std::queue<int> q;
  for(unsigned int i=0;i<topo.size();i++){
    if(dependencies(i)==0)// && flowdirs(i)!=NO_FLOW)  //Is it a peak?
      q.emplace(i);       
  }  //Yes.

  //Starting with the peaks, pass flow downstream
  progress.start(topo.size());
  while(!q.empty()){
    ++progress;

    const auto c = q.front();          //Copy focal cell from queue
    q.pop();                           //Clear focal cell from queue

    //Coordinates of downstream neighbour, if any
    const auto ndir = flowdirs(c); 

    int n = NO_FLOW;
    if(ndir!=NO_FLOW){  //TODO: Fix this monkey patching
      int x,y;
      topo.iToxy(c,x,y);
      const int nx = x+dx[ndir];
      const int ny = y+dy[ndir];
      n            = topo.xyToI(nx,ny);
      assert(n>=0);
    }

    //TODO: Might need this - could also check label
    //If downstream neighbour is the ocean, we drop our water into it and the
    //ocean is unaffected.
    // if(wtd(c)<=OCEAN_LEVEL){                                                    //I am confused how this is actually checking if downstream neighbour is the ocean. OCEAN_LEVEL = 0, so it looks more like here is a place to reset any accidentally negative wtd to 0? When coupled, negative wtd would actually be allowed, so we shouldn't do this
    //   wtd(c) = 0;
    //   continue;
    // }
  
    if (n == NO_FLOW){    //if this is a pit cell, move the water to the appropriate depression's water_vol.    
      if(wtd(c)>0){
        deps[label(c)].water_vol += wtd(c);
        wtd(c) = 0; //Clean up as we go
      }
    } else {                               //not a pit cell
      //If we have water, pass it downstream.
      if(wtd(c)>0){ //Groundwater can go negative, so it's important to make sure that we are only passing positive water around
        wtd(n) += wtd(c);  //Add water to downstream neighbour. This might result in filling up the groundwater table, which could have been negative
        wtd(c)  = 0;       //Clean up as we go
      }
  
      //Decrement the neighbour's dependencies. If there are no more dependencies,
      //we can process the neighbour.
      if(--dependencies(n)==0){                            //CHECK this is a hack! Something is wrong with the dependencies matrix, should never go below 0 but it sometimes does. 
        assert(dependencies(n)>=0);
        q.emplace(n);                   //Add neighbour to the queue
      }
    }
  }
  progress.stop();

  std::cerr<<"t FlowInDepressionHierarchy: Surface water = "<<timer.stop()<<" s"<<std::endl;
}



//When water overflows from one depression into another, this function ensures
//that chained overflows and infiltration take place.
//
//A depression has three places it can put the water it's received.
//The depression will try to stash water in these three places sequentially.
//  1. It can store the water in itself
//  2. It can overflow into its neighbouring depression (by following a geolink to that depression's leaf)
//  3. It can overflow into its parent
//
//Options (2) and (3) result in a recursive call. If there's enough water,
//eventually repeated calls to (3) will bring the function to the parent of the
//depression that originally called it (through its neighbour). At this point we
//stash the water in the parent and exit.
//
//Since we might end up calling the function again and again if there's a
//complex series of overflows, the `jump_table` argument holds the location of
//where the water ultimately ended up. Everything between the original node and
//this destination is then full which means that we only traverse each part of a
//filled hierarchy once.
//
//Note that since we only call this function on the leaf nodes of depressions
//the jump_table only needs to use leaves as its keys.
//
//@return The depression where the water ultimately ended up
template<class elev_t>
dh_label_t OverflowInto(
  const dh_label_t                            root,
  const dh_label_t                            stop_node,
  DepressionHierarchy<elev_t>                &deps,
  std::unordered_map<dh_label_t, dh_label_t> &jump_table,  //Shortcut from one depression to its next known empty neighbour
  double                                      extra_water
){
  auto &this_dep = deps.at(root);

  //TODO: Could simulate water running down flowpath into depression so that wtd
  //fills up more realistically

  if(root==OCEAN)                        //We've reached the ocean
    return OCEAN;                        //Time to stop: there's nowhere higher in the depression hierarchy

  //FIRST PLACE TO STASH WATER: IN THIS DEPRESSION

  //We've gone around in a loop and found the original node's parent. That means
  //it's time to stop. (This may be the leaf node of another metadepression, the
  //ocean, or a standard node.)
  if(root==stop_node){                   //We've made a loop, so everything is full
    if(this_dep.parent==OCEAN)           //If our parent is the ocean
      return OCEAN;                      //Then the extra water just goes away
    else                                 //Otherwise
      this_dep.water_vol += extra_water; //This node, the original node's parent, gets the extra water
    return stop_node;
  }

  if(this_dep.water_vol<this_dep.dep_vol){                                              //Can this depression hold any water?
    const double capacity = this_dep.dep_vol - this_dep.water_vol;                      //Yes. How much can it hold?
    if(extra_water<capacity){                                                           //Is it enough to hold all the extra water?
      this_dep.water_vol  = std::min(this_dep.water_vol+extra_water,this_dep.dep_vol);  //Yup. But let's be careful about floating-point stuff
      extra_water         = 0;                                                          //No more extra water
    } else {                                                                            //It wasn't enough to hold all the water
      this_dep.water_vol = this_dep.dep_vol;                                            //So we fill it all the way.
      extra_water       -= capacity;                                                    //And have that much less extra water to worry about
    }
  }

  if(extra_water==0)                                                                    //If there's no more extra water
    return root;                                                                        //Call it quits

  //TODO: Use jump table

  //Okay, so there's extra water and we can't fit it into this depression

  //SECOND PLACE TO STASH WATER: IN THIS DEPRESSION'S NEIGHBOUR
  //Maybe we can fit it into this depression's overflow depression!

  auto &pdep = deps.at(this_dep.parent);
  if(this_dep.odep==NO_VALUE){      //Does the depression even have such a neighbour? 
    if(this_dep.parent!=OCEAN && pdep.water_vol==0) //At this point we're full and heading to our parent, so it needs to know that it contains our water
      pdep.water_vol += this_dep.water_vol;
    return jump_table[root] = OverflowInto(this_dep.parent, stop_node, deps, jump_table, extra_water);  //Nope. Pass the water to the parent
  }

  //Can overflow depression hold more water?
  auto &odep = deps.at(this_dep.odep);
  if(odep.water_vol<odep.dep_vol){  //Yes. Move the water geographically into that depression's leaf.
    if(this_dep.parent!=OCEAN && pdep.water_vol==0 && odep.water_vol+extra_water>odep.dep_vol) //It might take a while, but our neighbour will overflow, so our parent needs to know about our water volumes
      pdep.water_vol += this_dep.water_vol + odep.dep_vol;           //Neighbour's water_vol will equal its dep_vol
    return jump_table[root] = OverflowInto(this_dep.geolink, stop_node, deps, jump_table, extra_water);
  }

  //Okay, so the extra water didn't fit into this depression or its overflow
  //depression. That means we pass it to this depression's parent.

  //If we've got here we have a neighbour, but we couldn't stash water in the
  //neighbour because it was full. So we need to see if our parent knows about
  //us.
  if(this_dep.parent!=OCEAN && pdep.water_vol==0)
    pdep.water_vol += this_dep.water_vol + odep.water_vol;

  //THIRD PLACE TO STASH WATER: IN THIS DEPRESSION'S PARENT
  return jump_table[root] = OverflowInto(this_dep.parent, stop_node, deps, jump_table, extra_water);
}



template<class elev_t>
void MoveWaterInDepHier(
  int                                         current_depression,
  DepressionHierarchy<elev_t>                &deps,
  std::unordered_map<dh_label_t, dh_label_t> &jump_table
){
  if(current_depression==NO_VALUE)
    return;

  auto &this_dep = deps.at(current_depression);

  //Visit child depressions. When these both overflow, then we spread water
  //across them by spreading water across their common metadepression
  MoveWaterInDepHier(this_dep.lchild, deps, jump_table);
  MoveWaterInDepHier(this_dep.rchild, deps, jump_table);

  //Catch depressions that link to the ocean through this one. These are special
  //cases because we will never spread water across the union of these
  //depressions and the current depressions: they only flow into the current
  //depression
  for(const auto c: this_dep.ocean_linked)
    MoveWaterInDepHier(c, deps, jump_table);

  //If the current depression is the ocean then at this point we've visited all
  //of its ocean-linked depressions (the ocean has no children). Since we do not
  //otherwise want to modify the ocean we quit here.
  if(current_depression==OCEAN)
    return;

  {
    const int lchild = this_dep.lchild;
    const int rchild = this_dep.rchild;

    //Only if both children are full should their water make its way to this
    //parent
    if(lchild!=NO_VALUE
      && deps.at(lchild).water_vol==deps.at(lchild).dep_vol
      && deps.at(rchild).water_vol==deps.at(rchild).dep_vol
    )
    this_dep.water_vol += deps.at(lchild).water_vol + deps.at(rchild).water_vol;
  }

  //Each depression has an associated dep_vol. This is the TOTAL volume of the
  //meta-depression including all of its children. This property answers the
  //question, "How much water can this meta-depression hold?"

  //Each depression also has an asoociated water volume called `water_vol`. This
  //stores the TOTAL volume of the water in the depression. That is, for each
  //depression this indicates how much water is in this depression and its
  //children. However, if a depression has sufficient volume to contain its
  //water, then its water_vol will not propagate up to its parent. In this way
  //we can distinguish between depressions whose water needs to be spread versus
  //metadepressions whose children might need water spread, but which will not
  //receive any spreading themselves.

  //We are overflowing the depression
  if(this_dep.odep == OCEAN){
    //If a depression overflows directly into an ocean then its odep is the
    //ocean and so is its parent.

    //The current depression's outlet is into the ocean. Since the ocean can
    //absorb an infinite amount of water without changing its water volume, we
    //simply set the amount of water contained in the current depression to be
    //either its water_vol (the depression doesn't overflow) or equal to its
    //depression volume (the excess water is thrown into the ocean, which is
    //unaffected).
    this_dep.water_vol = std::min(this_dep.water_vol,this_dep.dep_vol);
  } else if(this_dep.water_vol>this_dep.dep_vol) {
    //The neighbouring depression is not the ocean and this depression is
    //overflowing (therefore, all of its children are full)
    assert(this_dep.lchild==NO_VALUE || deps.at(this_dep.lchild).water_vol==deps.at(this_dep.lchild).dep_vol);
    assert(this_dep.rchild==NO_VALUE || deps.at(this_dep.rchild).water_vol==deps.at(this_dep.rchild).dep_vol);

    //The marginal volume of this depression is larger than what it can hold, so
    //we determine the amount that overflows, the "extra water".
    double extra_water = this_dep.water_vol - this_dep.dep_vol;
  
    //Now that we've figured out how much excess water there is, we fill this
    //depression to its brim. Note that we don't use addition or subtraction
    //here to ensure floating-point equality.
    this_dep.water_vol = this_dep.dep_vol;

    //OverflowInto will initially send water to this depression's neighbour's
    //leaf depression via the geolink. If everything fills up, the water will
    //end up in this depression's parent. So at this point we don't have to
    //worry about the extra water here any more.
    OverflowInto(this_dep.geolink, this_dep.parent, deps, jump_table, extra_water);

    assert(
         this_dep.water_vol==0 
      || this_dep.water_vol<=this_dep.dep_vol
      || (this_dep.lchild==NO_VALUE && this_dep.rchild==NO_VALUE) 
      || (
              this_dep.lchild!=NO_VALUE && this_dep.rchild!=NO_VALUE
           && deps.at(this_dep.lchild).water_vol<this_dep.water_vol
           && deps.at(this_dep.rchild).water_vol<this_dep.water_vol
         )
    );
  }

  //All overflowing depressions should by now have overflowed all the way down
  //to the ocean. We must now spread the water in the depressions by setting
  //appropriate values for wtd.
}

 

//Simple data structure to hold information needed to spread water in a filled
//metadepression.
class SubtreeDepressionInfo {
 public:
  //One of the depressions at the bottom of the meta-depression. We use this to
  //identify a pit cell from which to start flooding.
  int   leaf_label = -1;          
  //The metadepression containing all of the children. This metadepression is
  //guaranteed to be large enough to hold all of the water of its children plus
  //whatever exists only in the metadepression itself. We use this to determine
  //the water and depression volumes.
  int   top_label = -1;
  //Here we keep track of which depressions are contained within the
  //metadepression. This allows us to limit the spreading function to cells
  //within the metadepression.
  std::unordered_set<int> my_labels;
};



template<class elev_t, class wtd_t>
void FillDepressions(
  //Identifies a meta-depression through which water should be spread, leaf node
  //from which the water should be spread, valid depressions across which water
  //can spread, and the amount of water to spread
  SubtreeDepressionInfo             &stdi,  
  double                             water_vol, //Amount of water to spread around this depression
  const DepressionHierarchy<elev_t> &deps,  //Depression hierarchy
  const rd::Array2D<float>          &topo,  //Topographic data for calculating marginal volumes as we attempt to spread water
  const rd::Array2D<dh_label_t>     &label, //2D array in which each cell is labeled with the leaf depression it belongs to
  rd::Array2D<wtd_t>                &wtd    //Water table depth: we transfer water into this
){
  //Nothing to do if we have no water
  if(water_vol==0)
    return;

  //Hashset stores the ids of the cells we've visited. We don't want to use a 2D
  //array because the size of the DEM as a whole could be massive and that's a
  //lot of memory to allocate/deallocate each time this function is called. We
  //arbitrarily reserve enough space in the hashset for a few thousand items.
  //This should be large than most depressions while still being small by the
  //computer's standards.
  std::unordered_set<int> visited(2048);
 
  //Priority queue that sorts cells by lowest elevation first. If two cells are
  //of equal elevation the one added most recently is popped first. The ordering
  //of the cells processed by this priority queue need not match the ordering of
  //the cells processed by the depression hierarchy.
  rd::GridCellZk_high_pq<elev_t> flood_q;                          

  { //Scope to limit pit_cell
    //Cell from which we can begin flooding the meta-depression. Which one we
    //choose is arbitrary, since we will fill all of the leaf depressions and
    //intermediate metadepressions until and including when we reach the
    //depression identified by stdi.top_label
    const auto pit_cell    = deps.at(stdi.leaf_label).pit_cell;
    assert(pit_cell>=0);

    //We start flooding at the pit cell of the depression and work our way
    //upwards
    flood_q.emplace(
      pit_cell % topo.width(),
      pit_cell / topo.width(),
      topo(pit_cell)
    );

    visited.emplace(pit_cell);
  }

  //Cells whose wtd will be affected as we spread water around
  std::vector<int> cells_affected;

  //Stores the sum of the elevations of all of the cells in cells_affected. Used
  //for calculating the volume we've seen so far. (See explanation above or in
  //dephier.hpp TODO)
  double total_elevation = 0;

  while(!flood_q.empty()){
    const auto c = flood_q.top();
    flood_q.pop();

    //We keep track of the current volume of the depression by noting the total
    //elevation of the cells we've seen as well as the number of cells we've
    //seen.

    //TODO: Note that the current cell's above ground volume and wtd do not
    //contribute at all. This a choice that Kerry and Richard discussed. It is
    //as though there is a virtual water line coincident with the edge of the
    //current cell. No water infiltrates into this cell or is stored above it -
    //only cells previously visited are considered when doing volume
    //calculations.

    //Current volume of this subset of the metadepression. Since we might climb
    //over a saddle point, this value can occasionally be negative. It will be
    //positive by the time we need to spread the water around.
    const double current_volume = cells_affected.size()*topo(c.x,c.y) - total_elevation;

    //NOTE: If this is false by a small margin, then it's a floating point issue
    //and this should be adjusted to be >=-1e-6 and water_vol should be made 0
    assert(water_vol>=0); 

    //All the cells within this depression should have water table depths less
    //than or equal to zero because we have moved all of their water down slope
    //into the pit cell. Since we may already have filled other depressions
    //their cells are allowed to have wtd>0. Thus, we raise a warning if we are
    //looking at a cell in this unfilled depression with wtd>0.
    if(stdi.my_labels.count(label(c.x,c.y))==1 && wtd(c.x,c.y)>0)
      throw std::runtime_error("A cell was discovered in an unfilled depression with wtd>0!");

    //There are two possibilities:
    //1. The virtual water level exceeds slightly the height of the cell. The cell's water table then fills up as much as it can.
    //   The water surface is then level with the height of the cell.
    //
    //2. The water surface is below the height of the cell because there is sufficient topographic volume to hold all the water.
    //   In this case, the cell's water table is left unaffected.

    if(water_vol<=current_volume-wtd(c.x,c.y)){
      //The current scope of the depression plus the water storage capacity of
      //this cell is sufficient to store all of the water. We'll stop adding
      //cells and fill things now.

      //We will fill the depression so that the surface of the water is at this
      //elevation.
      double water_level;

      if(current_volume<water_vol){ //TODO: Check stdi.my_labels.count(label(c.x,c.y))==0 ?
        //The volume of water exceeds what we can hold above ground, so we will
        //stash as much as we can in this cell's water table. This is okay
        //because the above ground volume plus this cell's water table IS enough
        //volume (per the if-clause above).

        //Fill in as much of this cell's water table as we can
        const double fill_amount = water_vol - current_volume;
        assert(fill_amount>=0);
        wtd(c.x,c.y)   += fill_amount;
        water_vol -= fill_amount;   //Doesn't matter because we don't use water_vol anymore
        water_level     = topo(c.x,c.y);
      } else if (current_volume==water_vol){
        //The volume of water is exactly equal to the above ground volume so we
        //set the water level equal to this cell's elevation
        water_level = topo(c.x,c.y);
      } else {
        //The water volume is less than this cell's elevation, so we calculate
        //what the water level should be.

        //We have that Volume = (Water Level)*(Cell Count)-(Total Elevation)
        //rearranging this gives us:
        water_level = (water_vol+total_elevation)/cells_affected.size();
      }


      //Water level must be higher than (or equal to) the previous cell we looked at, but lower than (or equal to) the current cell
      assert(cells_affected.size()==0 || topo(cells_affected.back())<=water_level+FP_ERROR); 
      assert(topo(c.x,c.y)-water_level>=-1e-3);

      for(const auto c: cells_affected){
        assert(wtd(c)>=0);               //This should be true since we have been filling wtds as we go.
        if(water_level<topo(c)){
          assert(water_level>=topo(c)-FP_ERROR);
        }
        wtd(c) = water_level - topo(c);  //only change the wtd if it is an increase, here. We can't take water away from cells that already have it (ie reduce groundwater in saddle cells within a metadepression.)
        if(-FP_ERROR<=wtd(c) && wtd(c)<0)
          wtd(c) = 0;
        assert(wtd(c)>=0);
      }

      //We've spread the water, so we're done        
      return;
      
    }  else {
      //We haven't found enough volume for the water yet.

      //During the adding of neighbours neighbours might get added that are
      //lower than we are and belong to a different depression (notably, this
      //happens at the edge of a flat abuting an ocean). These cells will then
      //be popped and could be processed inappropriately. To prevent this, we
      //skip them here.
      if(stdi.my_labels.count(label(c.x,c.y))==0)  //CHECK. This was preventing cells that flowed to the ocean from allowing my depression volume to update. Is this way ok? Is this even needed?
        continue;

      //Okay, we're allow to add this cell's neighbours since this cell is part
      //of the metadepression.

      //Add this cell to those affected so that its volume is available for
      //filling.
      cells_affected.emplace_back(topo.xyToI(c.x,c.y));

      //Fill in cells' water tables as we go
      assert(wtd(c.x,c.y)<=0);
      water_vol += wtd(c.x,c.y);  //We use += because wtd is less than or equal to zero
      wtd(c.x,c.y)    = 0;             //Now we are sure that wtd is 0, since we've just filled it
      
      //Add the current cell's information to the running total
      total_elevation += topo(c.x,c.y);

      for(int n=1;n<=neighbours;n++){
        const int nx = c.x + dx[n]; //TODO ModFloor(x+dx[n],topo.width()); //Get neighbour's x-coordinate using an offset and wrapping
        const int ny = c.y + dy[n];                     //Get neighbour's y-coordinate using an offset
        if(!topo.inGrid(nx,ny))                         //Is this cell in the grid?
          continue;                                     //Nope: out of bounds.
        const int ni = topo.xyToI(nx,ny);               //Get neighbour's flat index
     
        //Ocean cells may be found at the edge of a depression. They might get
        //added to this list even if there are other, lower, cells within the
        //depression which have not yet been explored. This happens when a flat
        //abutts an ocean. The side of the flat near the ocean will see the
        //ocean and try to add it. The ocean would then be called instead of
        //more cells within the depression. Therefore, we do not add ocean
        //cells.

        //We must use the ocean level rather than the ocean label, or we will
        //mistakenly miss adding higher cells which belong to the ocean's depression
        //e.g. an escarpment before the ocean. 

        if(visited.count(ni)==0 && (label(nx,ny)!=OCEAN || topo(nx,ny)>OCEAN_LEVEL)){
          flood_q.emplace(nx,ny,topo(nx,ny));
          visited.emplace(ni);
        }
      }
    }
  }

  //Since we're in this function we are supposed to be guaranteed to be able to
  //fill our depression, since we only enter this function once that is true.
  //Therefore, if we've reached this point, something has gone horribly wrong
  //somewhere. :-(

  throw std::runtime_error("PQ loop exited without filling a depression!");
}



template<class elev_t>
SubtreeDepressionInfo FindDepressionsToFill(
  const int                          current_depression,    //Depression we are currently in
  const DepressionHierarchy<elev_t> &deps,                  //Depression hierarchy
  const rd::Array2D<float>          &topo,                  //Topographic data (used for determinining volumes as we're spreading stuff)
  const rd::Array2D<dh_label_t>     &label,                 //Array indicating which leaf depressions each cell belongs to
  rd::Array2D<float>                &wtd                    //Water table depth
){
  //Stop when we reach one level below the leaves
  if(current_depression==NO_VALUE)
    return SubtreeDepressionInfo();

  const auto& this_dep = deps.at(current_depression);

  //We start by visiting all of the ocean-linked depressions. They don't need to
  //pass us anything because their water has already been transferred to this
  //metadepression tree by MoveWaterInDepHier(). Similar, it doesn't mater what their leaf
  //labels are since we will never spread water into them.
  for(const auto c: this_dep.ocean_linked)
    FindDepressionsToFill(c, deps, topo, label, wtd);

  //At this point we've visited all of the ocean-linked depressions. Since all
  //depressions link to the ocean and the ocean has no children, this means we
  //have visited all the depressions and spread their water. Since we don't wish
  //to modify the ocean, we are done.
  if(current_depression==OCEAN)
    return SubtreeDepressionInfo();

  //We visit both of the children. We need to keep track of info from these
  //because we may spread water across them.
  SubtreeDepressionInfo left_info  = FindDepressionsToFill(this_dep.lchild, deps, topo, label, wtd);
  SubtreeDepressionInfo right_info = FindDepressionsToFill(this_dep.rchild, deps, topo, label, wtd);   

  SubtreeDepressionInfo combined;
  combined.my_labels.emplace(current_depression);
  combined.my_labels.merge(left_info.my_labels);
  combined.my_labels.merge(right_info.my_labels);

  combined.leaf_label = left_info.leaf_label;  //Choose left because right is not guaranteed to exist
  if(combined.leaf_label==NO_VALUE)            //If there's no label, then there was no child
    combined.leaf_label = current_depression;  //Therefore, this is a leaf depression

  combined.top_label = current_depression;

  //The water volume should never be greater than the depression volume because
  //otherwise we would have overflowed the water into the neighbouring
  //depression and moved the excess to the parent.
  assert(this_dep.water_vol<=this_dep.dep_vol);

  //Since depressions store their marginal water volumes, if a parent depression
  //has 0 marginal water volume, then both of its children have sufficient
  //depression volume to store all of their water. However, if our parent is an
  //ocean-link then we are guaranteed to be able to fill now because excess
  //water will have been transferred into the parent and we don't want to pool
  //the parent's water with our own (it might be at the bottom of a cliff).

  if(this_dep.water_vol<this_dep.dep_vol || this_dep.ocean_parent){
    assert(this_dep.water_vol<=this_dep.dep_vol);

    //If both of a depression's children have already spread their water, we do not
    //want to attempt to do so again in an empty parent depression. 
    //We check to see if both children have finished spreading water. 

    FillDepressions(combined, this_dep.water_vol, deps, topo, label, wtd);

    //At this point there should be no more water all the way up the tree until
    //we pass through an ocean link, so we pass this up as a kind of null value.
    return SubtreeDepressionInfo();
  } else {
    return combined;
  }
}



template<class elev_t, class wtd_t>
void FlowInDepressionHierarchy(
  const rd::Array2D<elev_t>     &topo,
  const rd::Array2D<dh_label_t> &label,
  const rd::Array2D<flowdir_t>  &flowdirs,
  DepressionHierarchy<elev_t>   &deps,
  rd::Array2D<wtd_t>            &wtd
){
  rd::Timer timer_overall;
  timer_overall.start();
  
  MoveWaterIntoPits(topo, wtd, label, deps, flowdirs);

  { 
    //Scope to limit `timer_overflow` and `jump_table`. Also ensures
    //`jump_table` frees its memory
    rd::Timer timer_overflow;
    timer_overflow.start();
    std::unordered_map<dh_label_t, dh_label_t> jump_table;
    MoveWaterInDepHier(OCEAN, deps, jump_table);
    std::cerr<<"t FlowInDepressionHierarchy: Overflow time = "<<timer_overflow.stop()<<std::endl;
  }

  //Sanity checks
  for(int d=1;d<(int)deps.size();d++){
    const auto &dep = deps.at(d);
    assert(dep.water_vol==0 || dep.water_vol<=dep.dep_vol);
    assert(dep.water_vol==0 || (dep.lchild==NO_VALUE && dep.rchild==NO_VALUE) || (dep.lchild!=NO_VALUE && deps.at(dep.lchild).water_vol<dep.water_vol));
    assert(dep.water_vol==0 || (dep.lchild==NO_VALUE && dep.rchild==NO_VALUE) || (dep.rchild!=NO_VALUE && deps.at(dep.rchild).water_vol<dep.water_vol));
  }

  std::cerr<<"p Finding filled..."<<std::endl;
  rd::Timer timer_filled;
  timer_filled.start();
  //This should check everything that is an immediate child of the ocean, so
  //we're supposed to hit all the depressions like this.
  FindDepressionsToFill(OCEAN,deps,topo,label,wtd);                              
  std::cerr<<"t FlowInDepressionHierarchy: Fill time = "<<timer_filled.stop()<<" s"<<std::endl;


  //TODO
  std::cerr<<"m Checking against master..."<<std::endl;
  rd::Array2D<float> wtd_master("wtd_master.dat", true);   //Recharge (Percipitation minus Evapotranspiration)
  for(unsigned int i=0;i<wtd.size();i++)
    assert(wtd(i)==wtd_master(i));
  std::cerr<<"m wtd field matches master!"<<std::endl;

  std::cerr<<"t FlowInDepressionHierarchy = "<<timer_overall.stop()<<" s"<<std::endl;
}

}

}

#endif
