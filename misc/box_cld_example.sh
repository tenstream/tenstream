#!/bin/bash
set -euo pipefail

make -j10 ex_pprts_box_cld
rm -rf edir.h5 ediff.h5 abso.h5
bin/ex_pprts_box_cld -show_edir hdf5:edir.h5 -show_ediff hdf5:ediff.h5 -show_abso hdf5:abso.h5 -solver 3_10 -pprts_solver_view -pprts_open_bc $@

ipython --pylab -c \
  'import h5py as H; \
  edir = H.File("edir.h5","r")["edir0"][:]; \
  ediff = H.File("ediff.h5","r")["ediff0"][:]; \
  abso = H.File("abso.h5","r")["abso0"][:]; \
  Nx,Ny = int(edir.shape[1]/2-.5), int(edir.shape[0]/2-.5); \
  print(Nx, Ny);\
  figure(figsize=(14,8));\
  subplot(4,3, 1); imshow(edir [Ny,:,:,0].T,interpolation="nearest");                xlabel("Nx"); ylabel("Nz");cb=colorbar();\
  subplot(4,3, 2); imshow(edir [:,Nx,:,0].T,interpolation="nearest");                xlabel("Ny"); ylabel("Nz");cb=colorbar();\
  subplot(4,3, 3); imshow(edir [:,:,-1,0],  interpolation="nearest", origin="lower");xlabel("Nx"); ylabel("Ny");cb=colorbar();cb.set_label("edir");\
  subplot(4,3, 4); imshow(ediff[Ny,:,:,0].T,interpolation="nearest");                xlabel("Nx"); ylabel("Nz");cb=colorbar();\
  subplot(4,3, 5); imshow(ediff[:,Nx,:,0].T,interpolation="nearest");                xlabel("Ny"); ylabel("Nz");cb=colorbar();\
  subplot(4,3, 6); imshow(ediff[:,:,-1,0],  interpolation="nearest", origin="lower");xlabel("Nx"); ylabel("Ny");cb=colorbar();cb.set_label("eup");\
  subplot(4,3, 7); imshow(ediff[Ny,:,:,1].T,interpolation="nearest");                xlabel("Nx"); ylabel("Nz");cb=colorbar();\
  subplot(4,3, 8); imshow(ediff[:,Nx,:,1].T,interpolation="nearest");                xlabel("Ny"); ylabel("Nz");cb=colorbar();\
  subplot(4,3, 9); imshow(ediff[:,:,-1,1],  interpolation="nearest", origin="lower");xlabel("Nx"); ylabel("Ny");cb=colorbar();cb.set_label("edn");\
  subplot(4,3,10); imshow(abso [Ny,:,:].T,interpolation="nearest");                xlabel("Nx"); ylabel("Nz");cb=colorbar();\
  subplot(4,3,11); imshow(abso [:,Nx,:].T,interpolation="nearest");                xlabel("Ny"); ylabel("Nz");cb=colorbar();\
  subplot(4,3,12); imshow(abso [:,:,-1],  interpolation="nearest", origin="lower");xlabel("Nx"); ylabel("Ny");cb=colorbar();cb.set_label("abso");\
  savefig("box_cld.pdf")'
