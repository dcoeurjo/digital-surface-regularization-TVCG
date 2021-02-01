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

//UI variables
float p_alpha  = 0.001;
float p_beta   = 1.0;
float p_gamma  = 0.1;
bool  p_clamp  = false;
float p_radius = 3.0;
int   p_nbSteps  = 50;
float p_dt       = 0.5;
float p_epsilon  = 0.00001;

//Polysocpe
polyscope::SurfaceMesh* primalSurf;
CountedPtr<SH3::BinaryImage> binary_image;
CountedPtr< SH3::DigitalSurface > surface;
SH3::SurfelRange surfels;
SH3::RealVectors normalsII;
SH3::KSpace K;

std::vector<std::vector<unsigned int>> faces;

template<typename Surf,typename Normals>
SHG3::RealVectors regularize(const Surf &surface,
                             const Normals &normals,
                             unsigned int nbsteps,
                             double alpha,
                             double beta,
                             double gamma,
                             bool clamp,
                             double dt,
                             double epsilon)
{
  DigitalSurfaceRegularization<SH3::DigitalSurface> regul(surface);
  regul.init(alpha,beta,gamma);
  auto surfelIndex = regul.getSurfelIndex();
  regul.attachNormalVectors([&](SH3::SCell &c){ return normals[ surfelIndex[c] ];});
  regul.enableVerbose();
  
  if (clamp)
    regul.regularize(nbsteps,dt,epsilon, DigitalSurfaceRegularization<SH3::DigitalSurface>::clampedAdvection);
  else
    regul.regularize(nbsteps,dt,epsilon);
  
  auto reg_pos = regul.getRegularizedPositions();
  auto posIndex = regul.getCellIndex();
  std::vector<RealPoint> newposindex;
  for(auto &surfel: surfels)
  {
    auto verts = SH3::getPrimalVertices(K, surfel, false );
    for(auto &v: verts)
      newposindex.push_back(reg_pos[ posIndex[v] ] );
  }
  return newposindex;
}

void doWorkNormals()
{
  //Computing some differential quantities
  auto params = SH3::defaultParameters() | SHG3::defaultParameters() |  SHG3::parametersGeometryEstimation();
  
  params("r-radius", p_radius);
  std::cout<<"Computing normals... "<<std::endl;
  normalsII = SHG3::getIINormalVectors(binary_image, surfels, params);
  std::cout<< "done."<<std::endl;
  primalSurf->addFaceVectorQuantity("II normal vectors", normalsII);
}

void doWork()
{
  if (normalsII.size()==0)
    doWorkNormals();
  SH3::Cell2Index c2i;
  std::cout<<"Computing regularization... "<<std::endl;
  auto newpos = regularize(surface, normalsII, p_nbSteps, p_alpha, p_beta, p_gamma, p_clamp, p_dt,p_epsilon);
  std::cout<< "done."<<std::endl;
 
  polyscope::registerSurfaceMesh("Regularized surface",newpos,faces);
}


void myCallback()
{
  ImGui::SliderFloat("Normal vector estimation radius (Integral Invariants)", &p_radius, 0.0, 5.0);
  if (ImGui::Button("Update normal vectors"))
    doWorkNormals();
  ImGui::SliderFloat("Data attachment term", &p_alpha, 0.00001, 1.0);
  ImGui::SliderFloat("Alignment term", &p_beta, 0.00001, 1.0);
  ImGui::SliderFloat("Fairness term", &p_gamma, 0.00001, 1.0);
  ImGui::SliderInt("Grad.: Number of steps", &p_nbSteps, 0, 100);
  ImGui::SliderFloat("Grad.: Initial learning rate", &p_dt, 0.1, 1.0);
  ImGui::SliderFloat("Grad.: Epsilon stopping criterion", &p_epsilon, 0, 0.1);
  ImGui::Checkbox("Grad.: Clamping", &p_clamp);
  
  if (ImGui::Button("Update"))
    doWork();
}


int main(int argc, char **argv)
{
  CLI::App app{"displayPTS"};
  std::string filename;
  app.add_option("-i,--input,1", filename, "Input VOL file")->required()->check(CLI::ExistingFile);
  CLI11_PARSE(app,argc,argv);
  
  polyscope::init();
  auto params = SH3::defaultParameters() | SHG3::defaultParameters() |  SHG3::parametersGeometryEstimation();
  
  auto h=1.; //gridstep
  params("surfaceComponents", "All");
  
  binary_image = SH3::makeBinaryImage(filename, params );
  K               = SH3::getKSpace( binary_image );
  surface         = SH3::makeDigitalSurface( binary_image, K, params );
  surfels         = SH3::getSurfelRange( surface, params );
  auto embedder   = SH3::getCellEmbedder( K );

  //Need to convert the faces
  faces.clear();
  std::vector<RealPoint> positions;
  unsigned int cpt=0;
  for(auto &surfel: surfels)
  {
    auto verts = SH3::getPrimalVertices(K, surfel, false );
    for(auto &v: verts)
      positions.push_back(embedder(v));
    
    std::vector<unsigned int> face={cpt, cpt+1, cpt+2,cpt+3};
    cpt+=4;
    faces.push_back(face);
  }
    
  primalSurf = polyscope::registerSurfaceMesh("Primal surface",positions, faces)->setEdgeWidth(1.0)->setEdgeColor({1.,1.,1.});

  polyscope::state::userCallback = myCallback;
  polyscope::show();
  return 0;
}
