/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systemes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2021/01/25
 */
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

#include <igl/grid.h>
#include <igl/dual_contouring.h>
#include <igl/marching_cubes.h>


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
bool  p_constrained = false;

//Polysocpe
polyscope::SurfaceMesh* primalSurf;
CountedPtr<SH3::BinaryImage> binary_image;
CountedPtr< SH3::DigitalSurface > surface;
CanonicSCellEmbedder<SH3::KSpace> sembedder;

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

void doWorkDualContouring()
{
  if (normalsII.size()==0)
    doWorkNormals();
  
  std::cout<<"Computing DC... "<<std::endl;
  Eigen::MatrixXd V;
  Eigen::MatrixXi Q,F;
  Eigen::MatrixXd mcV,mcN;
  Eigen::MatrixXi mcF;
  
  //Remap the normals
  std::map<Point, RealPoint> surfeltonormal;
  for(auto i=0;i < normalsII.size(); ++i)
  {
    auto p = sembedder(surfels[i]);
    Point surfcenter(2*p[0],2*p[1],2*p[2]);
    surfeltonormal[ surfcenter ] = normalsII[i];
  }
  // Grid parameters
  auto low = binary_image->domain().lowerBound();
  auto up  = binary_image->domain().upperBound();
  std::cout<<low<<" "<<up<<std::endl;
  const Eigen::RowVector3d min_corner(low[0],low[1],low[2]);
  const Eigen::RowVector3d max_corner(up[0],up[1],up[2]);
  const int s = 256;
  int nx = up[0] - low[0]+1;
  int ny = up[1] - low[1]+1;
  int nz = up[2] - low[2]+1;
  const Eigen::RowVector3d step =
  (max_corner-min_corner).array()/(Eigen::RowVector3d(nx,ny,nz).array()-1);
  // Sparse grid below assumes regular grid
  assert((step(0) == step(1))&&(step(0) == step(2)));
  
  // Dual contouring parameters
  bool constrained = p_constrainted;
  bool triangles = false;
  bool root_finding = false;
  
  //Following libIGL tutorial 715
  
  // Indicator function
  const auto & f_digital = [&](const Eigen::RowVector3d & x)
  {
    Point p  ((int)std::round(x[0]),
              (int)std::round(x[1]),
              (int)std::round(x[2]));
    if ((binary_image->domain().isInside(p)) && binary_image->operator()(p))
      return -1.0;
    else
      return 1.0;
  };
  
  // Simple finite difference gradients
  const auto & fd = [](
                       const std::function<double(const Eigen::RowVector3d&)> & ff,
                       const Eigen::RowVector3d & x)->Eigen::RowVector3d
  {
    const double eps = 1.0 ;
    Eigen::RowVector3d g;
    for(int c = 0;c<3;c++)
    {
      const Eigen::RowVector3d xp = x+eps*Eigen::RowVector3d(c==0,c==1,c==2);
      const double fp = ff(xp);
      const Eigen::RowVector3d xn = x-eps*Eigen::RowVector3d(c==0,c==1,c==2);
      const double fn = ff(xn);
      g(c) = (fp-fn)/(2*eps);
    }
    return g;
  };
  
  const auto  &f_grad_digital = [&fd,&f_digital](const Eigen::RowVector3d & x)->Eigen::RowVector3d
  {
    return fd(f_digital,x).normalized();
  };
  
  const auto  &f_grad_mapII = [&](const Eigen::RowVector3d & x)->Eigen::RowVector3d
  {
    Point surfcenter(2*x[0],2*x[1],2*x[2]);
    RealPoint normal = surfeltonormal[ surfcenter ];
    if (normal == RealPoint::zero)
    {
      std::cout << "Unknown surfel" << surfcenter << std::endl;
      exit(21);
    }
    Eigen::RowVector3d  n(normal[0], normal[1], normal[2]);
    return n;
  };
  
  igl::dual_contouring(f_digital,f_grad_mapII,min_corner,max_corner,nx,ny,nz,constrained,triangles,root_finding,V,Q);
  polyscope::registerSurfaceMesh("DC digital/II gradient",V,Q);
  
  igl::dual_contouring(f_digital,f_grad_digital,min_corner,max_corner,nx,ny,nz,constrained,triangles,root_finding,V,Q);
  polyscope::registerSurfaceMesh("DC digital/finite-diff gradient",V,Q);
  
  //Marching-Cubes
  Eigen::MatrixXd GV;
  igl::grid(Eigen::RowVector3i(nx,ny,nz),GV);
  Eigen::VectorXd Gf(GV.rows());
  igl::parallel_for(GV.rows(),[&](const int i)
                    {
    GV.row(i).array() *= (max_corner-min_corner).array();
    GV.row(i) += min_corner;
    Gf(i) = f_digital(GV.row(i));
  },1000ul);
  igl::marching_cubes(Gf,GV,nx,ny,nz,0,mcV,mcF);
  polyscope::registerSurfaceMesh("MC",mcV,mcF);
}

void myCallback()
{
  ImGui::SliderFloat("Normal vector estimation radius (Integral Invariants)", &p_radius, 0.0, 5.0);
  if (ImGui::Button("Compute normal vectors"))
    doWorkNormals();
  ImGui::Separator();
  ImGui::Text("Coefficient for each energy term:");
  ImGui::SliderFloat("Data attachment term", &p_alpha, 0.00001, 1.0);
  ImGui::SliderFloat("Alignment term", &p_beta, 0.00001, 1.0);
  ImGui::SliderFloat("Fairness term", &p_gamma, 0.00001, 1.0);
  ImGui::Separator();
  ImGui::Text("Gradient descent parameters:");
  ImGui::SliderFloat("Initial learning rate", &p_dt, 0.1, 1.0);
  ImGui::SliderInt("Number of steps", &p_nbSteps, 0, 100);
  ImGui::SliderFloat("Gradient norm stop", &p_epsilon, 0, 0.01);
  ImGui::Checkbox("Clamping", &p_clamp);
  
  if (ImGui::Button("Compute"))
    doWork();

  ImGui::Separator();
  ImGui::TextWrapped("For comparison, we also include the Dual Contouring and Marching Cubes approaches on the indicator function.");
  ImGui::TextWrapped("For DC, two functors for the gradient are given: one using finite difference of the indicator function, the other one uses a remapping of the Integral Invariant normal vectors.");
  
  ImGui::Checkbox("Force DC to clamp the vertices", &p_constrained);
  if (ImGui::Button("Compute DC/MC"))
    doWorkDualContouring();
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
