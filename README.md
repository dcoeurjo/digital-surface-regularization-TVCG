# Digital Surface Regularization with Guarantees

This repository is an implementation of the article *[Digital surface regularization with guarantees](https://perso.liris.cnrs.fr/david.coeurjolly/publication/dcoeurjotvcg21/)*,
David Coeurjolly, Jacques-Olivier Lachaud, Pierre Gueth, IEEE Transactionson Visualization and Computer Graphics, January 2021

# Building

The code is a graphical interface to the CPU regularization code available
in [DGtal](dgtal.org) using [polyscope](polyscope.run).

To compile the project, just clone (with its submodules) the repository:

```
git clone --recursive https://github.com/dcoeurjo/digital-surface-regularization-TVCG.git
```

and then use **cmake** to build the projects. For instance:

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release"
make
```


The tool can be used in commandline suing `regularization -i
bunny-64.vol` for instance.

Additional Input VOL files can be found in the
[VolGallery](https://github.com/dcoeurjo/VolGallery) repository.


# Citing


``` bibtex
@Article{dcoeurjoReg2021,
  author =       {David Coeurjolly, Jacques-Olivier Lachaud, Pierre
                  Gueth},
  title =        {Digital surface regularization with guarantees},
  journal =      {{IEEE} Transactions on Visualization and Computer Graphics},
  year =         2021,
  DOI =          {10.1109/tvcg.2021.3055242},
  note =         {to appear},
  mont =         jan,
  abstract =     {Voxel based modeling is a very attractive way to
                  represent complex multi-material objects. Beside
                  artistic choices of pixel/voxel arts, representing
                  objects as voxels allows efficient and dynamic
                  interactions with the scene. For geometry processing
                  purposes, many applications in material sciences,
                  medical imaging or numerical simulation rely on a
                  regular partitioning of the space with labeled
                  voxels. In this article, we consider a variational
                  approach to reconstruct interfaces in multi-labeled
                  digital images. This approach efficiently produces
                  piecewise smooth quadrangulated surfaces with some
                  theoretical stability guarantee. Non-manifold parts
                  at intersecting interfaces are handled naturally by
                  our model. We illustrate the strength of our tool
                  for digital surface regularization, as well as voxel
                  art regularization by transferring colorimetric
                  information to regularized quads and computing
                  isotropic geodesic on digital surfaces.}
                  }

```



