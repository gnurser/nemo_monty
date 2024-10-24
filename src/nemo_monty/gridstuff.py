from os.path import join as pjoin
import os, sys
import numpy as np
import numpy.ma as ma
import netCDF4
from argparse import ArgumentParser
from collections import namedtuple


# def make_2_slices(hbounds, wideslice=True):

#     if wideslice:
#         widehbounds = hbounds[:]
#         for n, l in enumerate(widehbounds):
#             if n % 2 == 0:  # lower limits
#                 widehbounds[n] -= 1
#             else:  # upper limits
#                 if l == -1:
#                     widehbounds[n] = None
#                 else:
#                     widehbounds[n] += 1
#         hboundslist = [hbounds, widehbounds]
#     else:
#         hboundslist = [hbounds]

#     def make_slice(bounds2d):
#         return (
#             slice(bounds2d[0], bounds2d[1], None),
#             slice(bounds2d[2], bounds2d[3], None),
#         )

#     return [make_slice(hb) for hb in hboundslist]


# def grid_from_fext(fext):
#     fextgrid = {
#         "t": {"T", "I", "P", "SIGI"},
#         "u": {"U", "U2"},
#         "v": {"V", "V2"},
#         "w": {"W", "W2", "LSPV"},
#         "f": {"PSI"},
#     }
#     for g, f in fextgrid.items():
#         if fext in f:
#             return g
#     sys.exit("no grid for fext %s found" % fext)


class DomainFiles:
    def __init__(self, domdir):
        self.domdir = domdir

    def __call__(self, fext):
        return pjoin(self.domdir, fext + ".nc")


class GridStuff:
    def __init__(self, grids):
        self.fextgrid = {"t": "T", "u": "U", "v": "V", "f": "PSI", "w": "W", "vw": "VW"}
        if grids is not None:
            self.grids = grids
        else:
            self.grids = self.fextgrid.keys()
        self.bdsgrid = {"t": "f", "u": "v", "v": "u", "f": "t", "w": "f", "vw": "u"}
        self.bdsgrid_disp = {
            "t": (0, 0),
            "u": (0, 1),
            "v": (1, 0),
            "f": (1, 1),
            "w": (0, 0),
            "vw": (1, 0),
        }
        self.wfile = open("gridinfo.txt", "w")

    @staticmethod
    def get_cdfvble(vble, grid):
        grid0 = {"u": "t", "v": "t", "f": "t", "vw": "w"}
        grid1 = {"w": "t", "vw": "v"}
        grid2 = {"vw": "w", "f": "u"}
        grid3 = {"vw": "w", "f": "t", "v": "t", "u": "t"}
        vbles = []
        if "mask" in vble:
            x = grid1.get(grid, grid) + vble
            vbles.append(x)
        elif "_0" in vble:
            x = vble[:-2] + grid0.get(grid, grid) + "_0"
            vbles.append(x)
            vbles.append(x.replace("0", "1d"))
        elif "3" in vble:
            x = vble + grid2.get(grid, grid)
            vbles.append(x)
            x = vble + grid0.get(grid, grid)
            vbles.append(x + "_0")
            vbles.append(x + "_1d")
        elif "gdep" in vble:
            x = vble + grid3.get(grid, grid)
            vbles.append(x)
            vbles.append(x + "_0")
            vbles.append(x + "_1d")
        else:
            x = vble + grid1.get(grid, grid)
            vbles.append(x)
        # print vbles
        return vbles

    def set_meshdict(self, run, meshvbles, grid, verbose=False):
        meshdict = {}
        meshdict["fext"] = self.fextgrid[grid]
        if isinstance(run, str):
            run = DomainFiles(run)
        for meshtype in ["mask", "mesh_hgr", "mesh_zgr"]:
            meshfile = run(fext=meshtype)
            if verbose:
                print("Reading from meshfile ", meshfile)
            else:
                print("Reading from meshfile ", meshfile, file=self.wfile)
            f = netCDF4.Dataset(meshfile)
            avail_vars = f.variables.keys()
            for vble in meshvbles:
                cdfvbles = self.get_cdfvble(vble, grid)
                if verbose:
                    print("Reading ", cdfvbles, vble)
                else:
                    print("Reading ", cdfvbles, vble, file=self.wfile)

                for cdfvble in cdfvbles:
                    if cdfvble in avail_vars:
                        meshdict[vble] = f.variables[cdfvble][...].squeeze()
                        break
            f.close()
        return meshdict

    def find_meshes(self, meshkeys, meshdir, verbose=False):
        meshkeys.append("fext")
        Mesh = namedtuple("Mesh", meshkeys)
        meshes = {}
        for grid in self.grids:
            meshdict = self.set_meshdict(meshdir, meshkeys, grid, verbose=verbose)
            try:
                meshes[grid] = Mesh(**meshdict)
            except:
                for attribute, value in meshdict.items():
                    print("{} : {}".format(attribute, value))
                sys.exit()
        self.wfile.close()
        print("full details of grids in file gridinfo.txt")
        return meshes

    def get_tra_grids(self, tracdict):
        self.tracdict = tracdict
        tracers = [x for y in tracdict.values() for x in y]
        tracers.sort()
        return tracers

    def find_grid(self, tracer):
        for k, v in self.tracdict.items():
            if tracer in v:
                return k
        else:
            sys.exit("Cannot find grid for tracer %s" % tracer)

    def find_fext(self, tracer):
        fexttrac = {}
        return fexttrac.get(tracer, self.fextgrid[self.find_grid(tracer)])


class RealGridStuff(GridStuff):
    @staticmethod
    def fitslice(Nddict, name, slice):
        Nd = Nddict[name]
        ndim = Nd.ndim
        ndim_real = 0
        for n in Nd.shape:
            if n > 1:
                ndim_real += 1
        ndslice = len(slice)
        if ndslice > ndim_real:
            sys.exit("real dim is %i for %s" % (ndim_real, name))
        else:
            slice = (
                (ndim - ndim_real) * (0,)
                + int(ndim_real > ndslice) * (Ellipsis,)
                + slice
            )
            return Nd[slice]

    def set_meshdict(self, run, meshvbles, grid, meshfile=None, verbose=False):
        hslice, wide_hslice = self.hslice, self.wide_hslice
        meshdict = {}
        meshdict["fext"] = self.fextgrid[grid]
        if isinstance(run, str):
            run = DomainFiles(run)
        if meshfile is None:
            meshfiles = ["mask", "mesh_hgr", "mesh_zgr"]
        else:
            meshfiles = [meshfile]
        for meshtype in meshfiles:
            meshpath = run(fext=meshtype)
            if meshpath[-6:] == ".nc.nc":
                meshpath = meshpath[:-3]

            if os.access(meshpath, os.R_OK):
                f = netCDF4.Dataset(meshpath)
                if verbose:
                    print("Reading from meshfile ", meshpath)
                else:
                    print("Reading from meshfile ", meshpath, file=self.wfile)
            else:
                sys.exit("no access to %s" % meshpath)

            avail_vars = f.variables.keys()
            print("In %s avail_vars are " % meshpath, avail_vars, file=self.wfile)
            if "area" in meshvbles:
                meshvbles.remove("area")
                meshvbles.append("area")
            for vble in meshvbles:
                if "bnds" in vble:
                    cdfvbles = self.get_cdfvble(vble[:-4], self.bdsgrid[grid])
                    for cdfvble in cdfvbles:
                        if cdfvble in avail_vars:
                            bds = self.fitslice(f.variables, cdfvble, wide_hslice)
                            print(vble, " has shape ", bds.shape, file=self.wfile)
                            nyp2, nxp2 = bds.shape
                            ny, nx = nyp2 - 2, nxp2 - 2
                            dtype = bds.dtype
                            meshdict[vble] = np.empty([ny, nx, 4], dtype=dtype)
                            dj, di = self.bdsgrid_disp[grid]
                            meshdict[vble][..., 0] = bds[dj : ny + dj, di : nx + di]
                            meshdict[vble][..., 1] = bds[
                                dj : ny + dj, di + 1 : nx + 1 + di
                            ]
                            meshdict[vble][..., 2] = bds[
                                dj + 1 : ny + 1 + dj, di + 1 : nx + 1 + di
                            ]
                            meshdict[vble][..., 3] = bds[
                                dj + 1 : ny + 1 + dj, di : nx + di
                            ]
                            if verbose:
                                print("Created mesh[ %s ].%s" % (grid, vble))
                            else:
                                print(
                                    "Created mesh[ %s ].%s" % (grid, vble),
                                    file=self.wfile,
                                )
                            break
                elif vble == "area" and meshtype == "mesh_hgr":
                    xdict = {}
                    for xvble in ["e1", "e2"]:
                        if xvble not in meshvbles:
                            cdfvble = self.get_cdfvble(xvble, grid)
                            varNd = f.variables[cdfvble]
                            xdict[xvble] = self.fitslice(f.variables, cdfvble, hslice)
                        else:
                            xdict[xvble] = meshdict[xvble]
                    meshdict[vble] = xdict["e1"] * xdict["e2"]
                    if verbose:
                        print("Created mesh[ %s ].%s" % (grid, vble))
                    else:
                        print("Created mesh[ %s ].%s" % (grid, vble), file=self.wfile)
                else:
                    cdfvbles = self.get_cdfvble(vble, grid)
                    for cdfvble in cdfvbles:
                        if cdfvble in avail_vars:
                            vbleNd = f.variables[cdfvble]
                            ndim_real = 0
                            for n in vbleNd.shape:
                                if n > 1:
                                    ndim_real += 1
                            if ndim_real == 1:
                                meshdict[vble] = f.variables[cdfvble][...].squeeze()[:]
                            else:
                                meshdict[vble] = self.fitslice(
                                    f.variables, cdfvble, hslice
                                )
                            print(
                                vble,
                                " has cdf variable ",
                                cdfvble,
                                " dtype ",
                                meshdict[vble].dtype,
                                " and shape ",
                                meshdict[vble].shape,
                                file=self.wfile,
                            )
                            if verbose:
                                print("Created mesh[ %s ].%s" % (grid, vble))
                            else:
                                print(
                                    "Created mesh[ %s ].%s" % (grid, vble),
                                    file=self.wfile,
                                )
                            break
            f.close()
        return meshdict

    def find_meshes(
        self,
        meshkeys,
        meshdir,
        hslice,
        wide_hslice=None,
        meshfile=None,
        meshbnds=False,
        verbose=False,
    ):
        meshkeys.append("fext")
        meshkeys = list(set(meshkeys))
        print(meshkeys)
        if meshbnds:
            if "glam" in meshkeys:
                meshkeys.append("glambnds")
            if "gphi" in meshkeys:
                meshkeys.append("gphibnds")
        Mesh = namedtuple("Mesh", meshkeys)
        meshes = {}
        self.hslice = hslice
        self.wide_hslice = wide_hslice
        for grid in self.grids:
            meshdict = self.set_meshdict(
                meshdir, meshkeys, grid, meshfile=meshfile, verbose=verbose
            )
            meshes[grid] = Mesh(**meshdict)
        self.wfile.close()
        print("full details of grids in file gridinfo.txt")
        return meshes
