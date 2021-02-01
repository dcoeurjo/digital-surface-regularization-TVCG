#include <iostream>
#include <fstream>
#include <string>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include <DGtal/geometry/surfaces/DigitalSurfaceRegularization.h>

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include "deps/CLI11/CLI11.hpp"


using namespace DGtal;
using namespace Z3i;

// Using standard 3D digital space.
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;

template<typename Surf,typename Norm>
SHG3::RealVectors regularize(const Surf &surface,
                             const Norm &normals,
                             const unsigned int nbsteps,
			     SH3::Cell2Index c2i )
{
  DigitalSurfaceRegularization<SH3::DigitalSurface> regul(surface);
  regul.init();
  auto surfelIndex = regul.getSurfelIndex();
  regul.attachNormalVectors([&](SH3::SCell &c){ return normals[ surfelIndex[c] ];});
  regul.regularize(nbsteps);

  auto reg_c2i = regul.getCellIndex();
  auto reg_pos = regul.getRegularizedPositions();
  auto ord_pos = reg_pos;
  for ( auto p : reg_c2i )
    {
      const auto pointel = p.first;
      const auto     idx = p.second;
      const auto new_idx = c2i[ pointel ];
      ord_pos[ new_idx ] = reg_pos[ idx ];
    }
  return ord_pos;
  //  return regul.getRegularizedPositions();
}
template<typename Surf,typename Norm>
SHG3::RealVectors regularizeClamp(const Surf &surface,
                             const Norm &normals,
                             const unsigned int nbsteps,
           SH3::Cell2Index c2i )
{
  DigitalSurfaceRegularization<SH3::DigitalSurface> regul(surface);
  regul.init();
  auto surfelIndex = regul.getSurfelIndex();
  regul.attachNormalVectors([&](SH3::SCell &c){ return normals[ surfelIndex[c] ];});
  regul.regularize(nbsteps,1.0,0.0001, DigitalSurfaceRegularization<SH3::DigitalSurface>::clampedAdvection);

  auto reg_c2i = regul.getCellIndex();
  auto reg_pos = regul.getRegularizedPositions();
  auto ord_pos = reg_pos;
  for ( auto p : reg_c2i )
    {
      const auto pointel = p.first;
      const auto     idx = p.second;
      const auto new_idx = c2i[ pointel ];
      ord_pos[ new_idx ] = reg_pos[ idx ];
    }
  return ord_pos;
  //  return regul.getRegularizedPositions();
}
template<typename Surf,typename Norm>
SHG3::RealVectors regularizeWOFairness(const Surf &surface,
                                       const Norm &normals,
                                       const unsigned int nbsteps,
           SH3::Cell2Index c2i )
{
  DigitalSurfaceRegularization<SH3::DigitalSurface> regul(surface);
  regul.init(0.001, 1.0,0.0);
  auto surfelIndex = regul.getSurfelIndex();
  regul.attachNormalVectors([&](SH3::SCell &c){ return normals[ surfelIndex[c] ];});
  regul.regularize(nbsteps);

  auto reg_c2i = regul.getCellIndex();
  auto reg_pos = regul.getRegularizedPositions();
  auto ord_pos = reg_pos;
  for ( auto p : reg_c2i )
    {
      const auto pointel = p.first;
      const auto     idx = p.second;
      const auto new_idx = c2i[ pointel ];
      ord_pos[ new_idx ] = reg_pos[ idx ];
    }
  return ord_pos;
}



int main()
{
  polyscope::init();
  
  auto params = SH3::defaultParameters() | SHG3::defaultParameters() |  SHG3::parametersGeometryEstimation();
  
  auto h=1.; //gridstep
  params( "polynomial", "goursat" )( "gridstep", h );//("minAABB",0);
  params("surfaceComponents", "All");

  /*auto implicit_shape  = SH3::makeImplicitShape3D  ( params );
  auto digitized_shape = SH3::makeDigitizedImplicitShape3D( implicit_shape, params );
  auto binary_image    = SH3::makeBinaryImage( digitized_shape, params );
  auto K               = SH3::getKSpace( params );
  */
  
  auto binary_image = SH3::makeBinaryImage("/Volumes/Stash/VolGallery/Stanford-bunny/bunny-128.vol", params );
  auto K               = SH3::getKSpace( binary_image );
    
  auto surface         = SH3::makeDigitalSurface( binary_image, K, params );
  auto embedder        = SH3::getCellEmbedder( K );
  auto surfels         = SH3::getSurfelRange( surface, params );
 
  
  //auto primalSurface   = SH3::makePrimalPolygonalSurface(c2i, surface);
  
  //Need to convert the faces
  std::vector<RealPoint> positions;
  unsigned int cpt=0;
  std::vector<std::vector<unsigned int>> facesbis;
  for(auto &surfel: surfels)
  {
    auto verts = SH3::getPrimalVertices(K, surfel, false );
    for(auto &v: verts)
      positions.push_back(embedder(v));
    
    std::vector<unsigned int> face={cpt, cpt+1, cpt+2,cpt+3};
    cpt+=4;
    facesbis.push_back(face);
  }
  
  SH3::Cell2Index c2i;
  
  auto digsurf = polyscope::registerSurfaceMesh("Primal surface",positions, facesbis);
  digsurf->setEdgeWidth(1.0);
  digsurf->setEdgeColor({1.,1.,1.});

  //Computing some differential quantities
  params("r-radius", 3.0);
  auto normalsII = SHG3::getIINormalVectors(binary_image, surfels, params);
  digsurf->addFaceVectorQuantity("II normal vectors", normalsII);
  
  //Reg
  DigitalSurfaceRegularization<SH3::DigitalSurface> regul( surface );
  regul.init();
  auto surfelIndex = regul.getSurfelIndex();
  regul.attachNormalVectors([&](SH3::SCell &c){ return normalsII[ surfelIndex[c] ];} );
  regul.regularize();
 
  auto  newpos = regul.getRegularizedPositions();
  
  auto posIndex = regul.getCellIndex();
  std::vector<RealPoint> newposindex;
  for(auto &surfel: surfels)
  {
    auto verts = SH3::getPrimalVertices(K, surfel, false );
    for(auto &v: verts)
      newposindex.push_back(newpos[ posIndex[v] ] );
  }
  
  auto reg = polyscope::registerSurfaceMesh("Regularized surface", newposindex, facesbis);
  reg->setEdgeWidth(1.0);
  reg->setEdgeColor({1.,1.,1.});
  
  
  //Reg
  DigitalSurfaceRegularization<SH3::DigitalSurface> regul2( surface );
  regul2.init();
  auto surfelIndex2 = regul2.getSurfelIndex();
  regul2.attachNormalVectors([&](SH3::SCell &c){ return normalsII[ surfelIndex2[c] ];} );
  regul2.regularize(200,1.0,0.0001, DigitalSurfaceRegularization<SH3::DigitalSurface>::clampedAdvection );
  
  auto  newpos2 = regul2.getRegularizedPositions();
  std::vector<RealPoint> newposindex2;
  for(auto &surfel: surfels)
  {
    auto verts = SH3::getPrimalVertices(K, surfel, false );
    for(auto &v: verts)
      newposindex2.push_back(newpos2[ posIndex[v] ] );
  }
  auto reg2 = polyscope::registerSurfaceMesh("Regularized surface (clamp)", newposindex2, facesbis);
  reg2->setEdgeWidth(1.0);
  reg2->setEdgeColor({1.,1.,1.});
  
  
  //Reg without fairness
  DigitalSurfaceRegularization<SH3::DigitalSurface> regulWOF( surface );
  regulWOF.init( 0.001, 1.0, 0.00000000000001);
  auto surfelIndexWOF = regulWOF.getSurfelIndex();
  regulWOF.attachNormalVectors([&](SH3::SCell &c){ return normalsII[ surfelIndex2[c] ];} );
  regulWOF.regularize(200,1.0,0.0001);
  
  auto  newposWOF = regulWOF.getRegularizedPositions();
  std::vector<RealPoint> newposindexWOF;
  for(auto &surfel: surfels)
  {
    auto verts = SH3::getPrimalVertices(K, surfel, false );
    for(auto &v: verts)
      newposindexWOF.push_back(newposWOF[ posIndex[v] ] );
  }
  auto regWOF = polyscope::registerSurfaceMesh("Regularized surface (without fairness)", newposindexWOF, facesbis);
  regWOF->setEdgeWidth(1.0);
  regWOF->setEdgeColor({1.,1.,1.});
  
  //Reg without fairness with clamp
  DigitalSurfaceRegularization<SH3::DigitalSurface> regulWOFClamp( surface );
  regulWOFClamp.init( 0.001, 1.0, 0.00000000000001);
  auto surfelIndexWOFClamp = regulWOFClamp.getSurfelIndex();
  regulWOFClamp.attachNormalVectors([&](SH3::SCell &c){ return normalsII[ surfelIndex2[c] ];} );
  regulWOFClamp.regularize(200,1.0,0.0001,DigitalSurfaceRegularization<SH3::DigitalSurface>::clampedAdvection );
  
  auto  newposWOFClamp = regulWOFClamp.getRegularizedPositions();
  std::vector<RealPoint> newposindexWOFClamp;
  for(auto &surfel: surfels)
  {
    auto verts = SH3::getPrimalVertices(K, surfel, false );
    for(auto &v: verts)
      newposindexWOFClamp.push_back(newposWOFClamp[ posIndex[v] ] );
  }
  auto regWOFClamp = polyscope::registerSurfaceMesh("Regularized surface (without fairness, with clamp)", newposindexWOFClamp, facesbis);
  regWOFClamp->setEdgeWidth(1.0);
  regWOFClamp->setEdgeColor({1.,1.,1.});
  
  
  //Reg without align
  DigitalSurfaceRegularization<SH3::DigitalSurface> regulWOA( surface );
  regulWOA.init( 0.001, 0.00000000000000001, 0.005);
  auto surfelIndexWOA = regulWOA.getSurfelIndex();
  regulWOA.attachNormalVectors([&](SH3::SCell &c){ return normalsII[ surfelIndex2[c] ];} );
  regulWOA.regularize(200,1.0,0.0001);
  
  auto  newposWOA = regulWOA.getRegularizedPositions();
  std::vector<RealPoint> newposindexWOA;
  for(auto &surfel: surfels)
  {
    auto verts = SH3::getPrimalVertices(K, surfel, false );
    for(auto &v: verts)
      newposindexWOA.push_back(newposWOA[ posIndex[v] ] );
  }
  auto regWOA = polyscope::registerSurfaceMesh("Regularized surface (without align)", newposindexWOA, facesbis);
  regWOA->setEdgeWidth(1.0);
  regWOA->setEdgeColor({1.,1.,1.});
  
  
  
  polyscope::show();
  return 0;
}
