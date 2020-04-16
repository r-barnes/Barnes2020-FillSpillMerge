#pragma once

#include <richdem/common/Array2D.hpp>
#include <richdem/common/timer.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <richdem/common/grid_cell.hpp>
#include <richdem/common/constants.hpp>

namespace richdem {

/**
  @brief  Fills all pits and removes all digital dams from a DEM, but faster
  @author Richard Barnes (rbarnes@umn.edu)

    Priority-Flood starts on the edges of the DEM and then works its way
    inwards using a priority queue to determine the lowest cell which has
    a path to the edge. The neighbours of this cell are added to the priority
    queue if they are higher. If they are lower, they are raised to the
    elevation of the cell adding them, thereby filling in pits. The neighbors
    are then added to a "pit" queue which is used to flood pits. Cells which
    are higher than a pit being filled are added to the priority queue. In this
    way, pits are filled without incurring the expense of the priority queue.

  @param[in,out]  &elevations   A grid of cell elevations

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM. Note that the _NoData_ value is assumed to
       be a negative number less than any actual data value.

  @post
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM.
    2. **elevations** contains no landscape depressions or digital dams.

  @correctness
    The correctness of this command is determined by inspection. (TODO)
*/
template <Topology topo, class elev_t>
void PriorityFlood_Barnes2014_OceanInit(
  Array2D<elev_t> &elevations,
  const elev_t ocean_level
){
  GridCellZ_pq<elev_t> open;
  std::queue<GridCellZ<elev_t> > pit;
  uint64_t processed_cells = 0;
  uint64_t pitc            = 0;
  ProgressBar progress;

  RDLOG_ALG_NAME << "Priority-Flood (Improved)";
  RDLOG_CITATION << "Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024";
  RDLOG_CONFIG   <<"topology = "<<TopologyName(topo);

  const int *const dx   = topo == Topology::D8?d8x:topo==Topology::D4?d4x:NULL;
  const int *const dy   = topo == Topology::D8?d8y:topo==Topology::D4?d4y:NULL;
  const int        nmax = topo == Topology::D8?  8:topo==Topology::D4?  4:   0;

  RDLOG_PROGRESS << "Setting up boolean flood array matrix...";
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);

  RDLOG_MEM_USE<<"Priority queue requires approx = "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
               <<"MB of RAM.";

  RDLOG_PROGRESS<<"Adding cells to the priority queue...";

  for(int y=0;y<elevations.height();y++)
  for(int x=0;x<elevations.height();x++){
    if(elevations.isNoData(x,y) || elevations(x,y)==ocean_level){
      open.emplace(x,y,elevations(x,y));
      closed(x,y)=true;
    }
  }

  RDLOG_PROGRESS<<"Performing the improved Priority-Flood...";
  progress.start( elevations.size() );
  while(open.size()>0 || pit.size()>0){
    GridCellZ<elev_t> c;
    if(pit.size()>0){
      c=pit.front();
      pit.pop();
    } else {
      c=open.top();
      open.pop();
    }
    processed_cells++;

    for(int n=1;n<=nmax;n++){
      int nx=c.x+dx[n];
      int ny=c.y+dy[n];
      if(!elevations.inGrid(nx,ny)) continue;
      if(closed(nx,ny))
        continue;

      closed(nx,ny)=true;
      if(elevations(nx,ny)<=c.z){
        if(elevations(nx,ny)<c.z){
          ++pitc;
          elevations(nx,ny)=c.z;
        }
        pit.push(GridCellZ<elev_t>(nx,ny,c.z));
      } else
        open.emplace(nx,ny,elevations(nx,ny));
    }
    progress.update(processed_cells);
  }
  RDLOG_TIME_USE<<"Succeeded in "<<std::fixed<<std::setprecision(1)<<progress.stop()<<" s";
  RDLOG_MISC    <<"Cells processed = "<<processed_cells;
  RDLOG_MISC    <<"Cells in pits = "  <<pitc;
}

}