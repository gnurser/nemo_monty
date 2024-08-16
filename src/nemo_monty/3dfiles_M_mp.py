#!/usr/bin/env python
# -*- coding: utf-8 -*-
from os.path import join as pjoin
import os
import sys
import numpy as np
import numpy.ma as ma
import netCDF4
from argparse import ArgumentParser
from glob import glob
import time
from numba import jit
import math

from nemo_monty.gridstuff import RealGridStuff
import nemo_monty.findnemo as findnemo
from nemo_eos.nemo_rho import eos
from nemo_monty.interp import interp
import bp

#@jit
def tracer_interpolate(tr,k_below_s,r_above_s,active,nopython=True):
    JM,IM = active.shape
    tr_s = np.zeros([JM,IM],dtype=tr.dtype)

    for j in range(JM):
        for i in range(IM):
            if active[j,i]:
                k = k_below_s[j,i]
                ra = r_above_s[j,i]
                tr_s[j,i] = ra*tr[k-1,j,i] + (1.-ra)*tr[k,j,i]

    return tr_s


class Data(object):
    def __init__(self,name):
        # if isinstance(name,basestring):
        if isinstance(name,str):
            self.name=name
        elif isinstance(name,netCDF4.Variable):
            self.name=name.name

class InArgs(object):
    def __init__(self,arguments,vdict={},start=[],end=[]):
        if not isinstance(start, (list, tuple, set)):
            start = [start]
        if not isinstance(end, (list, tuple, set)):
            end = [end]
        self.start = start
        self.end = end

        for x,y in vdict.items():
            if not isinstance(y, (list, tuple, set)):
                vdict[x] = {y}
            else:
                vdict[x] = set(y)

        arguments = set([self.strip(a) for a in arguments] + arguments)
        for a in arguments.copy():
            if a in vdict.keys():
                arguments |= vdict[a]
        arguments = set([self.strip(a) for a in arguments]) | arguments
        self.arguments = arguments

    def strip(self,vble):
        for s in self.start:
            # print s
            if vble.startswith(s):
                l = len(s)
                vble = vble[l:]
        for e in self.end:
            if vble.endswith(e):
                l = len(e)
                vble = vble[:-l]
        return vble

    def __call__(self,vname):
        if vname in self.arguments:
            return True
        else:
            return False

def data_like(d0,name):
    data = Data(name)
    d0D = d0.__dict__
    data._FillValue = None
    for att in set(['dimensions','coordinates','_FillValue','grid','online_operation']) & set(d0D.keys()):
        data.__setattr__(att,d0D[att])
    return data


def change_str(in_str,**kwargs):
    for old,new in kwargs.items():
        in_str = in_str.replace(old,new)

    return in_str


nemo_mean_names = {'T':'votemper','mld':'somxl010','ssh':'sossheig','S':'vosaline',
             'Hin':'sohefldo','EmP':'sowaflup','EmPs':'sowaflcd',
             'sss':'sosaline','sst':'sosstsst',
             'u':'vozocrtx','age':'Age','PWT':'PWT','lspv':'volspv','sigi':'vosigmai'}
nemo_restart_names = {'T':'tn','ssh':'ssh_m','S':'sn',
             'sss':'sss_m','sst':'sstm',
             'u':'un','hice':'hicif','hsnow':'hsnif','aice':'frld'}
nemo_names = {'T':{'votemper','tn','potemp','thetao_con'},'mld':{'somxl010','somxzint1'},'ssh':{'sossheig','ssh_m','ssh','zos'},
                  'S':{'vosaline','sn','salin','so_abs'},
             'Hin':{'sohefldo','hfds'},'EmP':{'sowaflup'},'EmPs':{'sowaflcd'},
             'sss':{'sosaline','sss_m', 'sss','sos_abs'},'sst':{'sosstsst','sstm','sst','tos_con'},
             'u':{'vozocrtx','un'},'age':{'Age'},'PWT':{'PWT'},
                  'lspv':{'volspv'},'sigi':{'vosigmai'},
                  'hice':{'hicif'},'hsnow':{'hsnif'},'aice':{'frld'} }
nemo_dimensions ={'T':3,'mld':2,'ssh':2,'S':3,
             'Hin':2,'THin':2,'EmP':2,'EmPs':2,
             'sss':2,'sst':2,
             'u':3,'age':3,'lspv':3,'sigi':3,'PWT':3}


class FextTrac(object):
    def __init__(self,vtype):
        if vtype=='mean':
            self.ftall = {'T':{'ssh','T','mld','S','Hin','EmP','EmPs','sst','sss'},
                            'I':{'hice','hsnow','aice'},
                            'F':{'barstf'},'U':{'u'},'V':{'v'},'P':{'age','PWT'},
                            'LSPV':{'lspv'},'SIG2':{'sigi'}}
        elif vtype=='restart':
            self.ftall = {'restart':{'ssh','T','mld','S','Hin','EmP','EmPs','sst','sss','u','v'},
                            'restart_ice':{'hice','hsnow','aice'},
                            'restart_trc':{'age','PWT'}}

        self.keytype = 'fext'

    def from_tracer(self,t):
        for g,allt in self.ftall.items():
            if t in allt or nemo_mean_names.get(t) in allt or  nemo_restart_names.get(t) in allt:
                return g
        sys.exit('no %s for tracer %s found' % (self.keytype,t))

    def get_tracdict(self,tracers):
        ftdict = {}
        if not tracers:
            return ftdict
        tracset = set(tracers)
        # print 'tracers are ',' '.join(tracset)
        # print '\n',self.ftall.keys()
        for f in self.ftall.keys():
            gtracers = {t for t in tracset if self.from_tracer(t)==f}
            if gtracers:
                ftdict[f] = gtracers
                tracset -= gtracers
        if tracset:
            sys.exit('no %s for tracers %s found' % (self.keytype,' '.join(tracset)))
        else:
            for k,v in ftdict.items():
                print('{%s : %s}' % (k,' '.join(v)), end='')
            #print('\n')
            return ftdict


class GridTrac(FextTrac):
    def __init__(self):
        self.ftall = {'t':{'ssh','T','mld','S','Hin','EmP','EmPs','sst','sss',
                           'hice','hsnow','aice','age','PWT','sigi'},
               'f':{'barstf'},'u':{'u'},'v':{'v'},'w':{'lspv'}}
        self.keytype = 'grid'


class DCDF4(object):
    @staticmethod
    def grid_from_depth(dimlist):
        depthlist = [x for x in dimlist if 'dep' in x]
        if depthlist:
            depth, = depthlist
            return depth[-1]
        else:
            return None

    @staticmethod
    def change_tup(in_tuple,**kwargs):
        in_list = list(in_tuple)
        for i,d in enumerate(in_list):
            for old,new in kwargs.items():
                if old in d:
                    in_list[i] = new

        return tuple(in_list)

    @classmethod
    def set_gridtrac(cls,gridtrac_instance):
        cls.gridtrac = gridtrac_instance

    @classmethod
    def set_default_slice(cls,default_slice):
        cls.default_slice = default_slice

    def __init__(self, f_or_fname, grid=None, slice=None, checkmask=False):
        if isinstance(f_or_fname,netCDF4.Dataset):
            f = f_or_fname
            self.f_keep = True
        elif isinstance(f_or_fname,str):
        # elif isinstance(f_or_fname,basestring):
            f = netCDF4.Dataset(f_or_fname)
            f.set_auto_mask(False)
            print('Opened %s' % f_or_fname)
            self.f_keep = False
        else:
            sys.exit('neither string nor netCDF4.Dataset')

        self.f = f
        self.dimensions = f.dimensions
        self.fv = f.variables
        self.slice = slice
        if grid is None:
            grid =  self.grid_from_depth(self.dimensions.keys())
        self.grid = grid
        self.checkmask = checkmask

    def __getitem__(self,slice):
        self.slice = slice

    def getNd(self,name):
        keys = self.fv.keys()
        if name not in keys:
            nemonames = nemo_names.get(name)
            if not nemonames:
                sys.exit('%s not in file or in nemonames' % name)
            else:
                nemokeys = nemonames & set(keys)
                if not nemokeys:
                    sys.exit('nemonames %s not in file' % ' '.join(nemonames))
                else:
                    nemoname, = nemokeys
        else:
            nemoname = name

        return self.fv[nemoname]

    def update(self,data):
        t0 = time.time()
        Nd = self.getNd(data.name)
        if self.slice is None:
            slice = self.default_slice
        else:
            slice = self.slice
        x = Nd[slice]

        # if data.name=='EmP':
        #     if  isinstance(x,ma.masked_array):
        #         x = x.data
        if hasattr(data, 'mask'):
                data.nos = ma.masked_where(data.mask,x)
        else:
            if not isinstance(x,ma.masked_array):
                if hasattr(Nd,'_FillValue'):
                    fillvalue = Nd._FillValue
                    data.nos = ma.masked_equal(x,fillvalue)
                else:
                    print('warning: data %s could not be masked from ._FillValue' % (data.name))
                    fillvalue = x.ravel()[-1]
                    print('trying masking by last value %g' % fillvalue)
                    data.nos = ma.masked_equal(x,fillvalue)
            else:
                data.nos = x
        t1 = time.time()
        print('reading %s' % data.name,'slice= ',slice,' shape= ',x.shape,
         ' took %f seconds' % (t1-t0))

    def get_tracer(self,name,grid=None,dtype=None,meshes=None):
        t0 = time.time()
        data = Data(name)
        Nd = self.getNd(name)
        NdD = Nd.__dict__
        data.dimensions = self.change_tup(Nd.dimensions,dep=u'z')
        # split so depth doesn't get changed into time-counter
        data.dimensions = self.change_tup(data.dimensions,t=u'time_counter')
        for g in (grid, NdD.get('grid'), self.gridtrac.from_tracer(name), self.grid_from_depth(data.dimensions),self.grid):
            if g is not None:
                data.grid = g
                break
        else:
            sys.exit('could not find grid for %s' % name)

        data._FillValue = None
        for att in set(['coordinates','units','long_name','_FillValue',
                        'online_operation','standard_name','short_name']) & set(NdD.keys()):
            data.__setattr__(att,NdD[att])

        if 'coordinates' in data.__dict__.keys():
            data.coordinates = change_str(data.coordinates,
                                    time_counter=u'years',
                                    nav_lon=u'glam%s' % data.grid,
                                    nav_lat=u'gphi%s' % data.grid,
                                    depth=u'gdep%s' % data.grid)
        elif 'x' in data.dimensions and 'y' in data.dimensions:
            data.coordinates = u'glam%s gphi%s' % (2*(data.grid,))

        if self.slice is None:
            slice = self.default_slice
        else:
            slice = self.slice
        data.nos = Nd[slice]

        data.otherdims = {}
        for i,d in enumerate(data.dimensions):
            if d not in ('time_counter','z','y','x'):
                data.otherdims[d] = self.getNd(d)[slice[i]]

        def get_mask( ):
            if data.nos.ndim == 2:
                mask = meshes[data.grid].maskutil
            elif data.nos.ndim == 3:
                mask = meshes[data.grid].mask

            return ~(mask.astype(bool))

        # if name == 'EmP':
        #     oldmask = data.nos.mask
        #     data.nos = data.nos.data
        # else:
        #     oldmask = None
        # print('checkmask=', self.checkmask)

        if not self.checkmask:
            if not isinstance(data.nos, ma.masked_array):
                dmax = data.nos.max()
                print('%s  ' % name,end='')
                if dmax > 1.e12:
                    data.nos = ma.masked_equal(data.nos,dmax)
                    print('dmax loop ',dmax,data.nos.max())
                    data.nos.mask += data.nos.data==0.
                elif data._FillValue is not None:
                    data.nos = ma.masked_equal(data.nos,data._FillValue)
                    print('FV loop ', data._FillValue,data.nos.max())
                else:
                    data.nos = ma.masked_equal(data.nos,0.)
                    print('zero loop ', data.nos.max())

                # if oldmask is not None:
                #     # data.nos[oldmask] = ma.masked
                #     data.nos.mask += oldmask

                data.mask = data.nos.mask

        elif meshes is not None:
            data.mask = get_mask()
            if not isinstance(data.nos, ma.masked_array):
                data.nos = ma.masked_where(data.mask,data.nos)
            else:
                data.nos = ma.masked_where(data.mask, data.nos.data)

        t1 = time.time()
        print('reading %s' % data.name,'slice= ',slice,' shape= ',data.nos.shape,
          ' took %f seconds' % (t1-t0))
        return data

    def __call__(self,tracers, meshes=None):
        P = {}
        for tracer in tracers:
            P[tracer] = self.get_tracer(tracer, meshes=meshes)
        return P


    def __del__(self):
        if not self.f_keep: self.f.close()


class Create3DCDF():
    def __init__(self,grids,meshes,outpath='scalars.nc',
                 time_index_name='years',dims=['t','z','y','x'], other_dims_dict = {},density=None):
        self.time_index_name = time_index_name
        f = netCDF4.Dataset(outpath,'w')
        self.f = f

        # ny,nx = meshes[grids[0]].glam.shape
        ny,nx = meshes['t'].glam.shape
        # ny -=2
        # nx -=2
        if 'y' in dims:
            f.createDimension('y',ny)
            yNd = f.createVariable('y',np.int32,('y',))
            yNd[:] = np.arange(ny) + 1
            yNd.standard_name = "j-index"
            yNd.long_name = "NEMO j-index (Fortran numbering)"

        if 'x' in dims:
            f.createDimension('x',nx)
            xNd = f.createVariable('x',np.int32,('x',))
            xNd[:] = np.arange(nx) + 1
            xNd.standard_name = "i-index"
            xNd.long_name = "NEMO i-index (Fortran numbering)"

        for d,dvals in other_dims_dict.items():
            nd = len(dvals)
            f.createDimension(d,nd)
            xNd = f.createVariable(d,np.float32,(d,))
            xNd[:] = np.array(dvals)
            xNd.standard_name = d
            xNd.long_name = d
            if d=='sigma':
                xNd.units = 'kg m^-3'
                xNd.positive = 'down'
                xNd._CoordinateAxisType = "Depth"

        # if 'x' in dims or 'y' in dims:
        #     if 'x' in dims and 'y' in dims:
        #         nv = 4
        #     else:
        #         nv = 2
        #     htuple = tuple([d for d in dims.split() if d in ('x', 'y')])
        #     f.createDimension('vertices', nv)

        htuple = tuple([d for d in list(dims) if d in ('y', 'x')])
        if 'x' in dims and 'y' in dims:
            nv = 4
            f.createDimension('vertices', nv)
            vtuple = ('vertices',)
        else:
            vtuple = ()

        if 'z' in dims:
            nz, = meshes[grids[0]].gdep_0.shape
            f.createDimension('z',nz)
            xNd = f.createVariable('z',np.int32,('z',))
            xNd[:] = np.arange(nz) + 1
            xNd.standard_name = "k-index"
            xNd.positive = 'down'
            xNd.long_name = "NEMO k-index (Fortran numbering)"

        #print(dims)

        if 't' in dims:
            f.createDimension('time_counter', None)
            xNd = f.createVariable('time_counter',np.int32,('time_counter',))
            xNd.standard_name = "time-index"
            xNd.long_name = "NEMO dataset number from start"
            self.TCNd = xNd
            self.TC = 0

            time_indexNd = f.createVariable(time_index_name,np.float32,('time_counter',))
            self.time_indexNd = time_indexNd
            if time_index_name == 'years':
                time_indexNd.units = 'common_year'
                time_indexNd.calendar = 'noleap'
                time_indexNd.standard_name = "years"
                time_indexNd.long_name = "years"

                self.set_time_index(time_index_name)

            dtime_indexNd = f.createVariable('time_interval',np.float64,('time_counter',))
            self.dtime_indexNd = dtime_indexNd
            if time_index_name == 'years':
                dtime_indexNd.units = 'days'
                dtime_indexNd.standard_name = "time interval"
                dtime_indexNd.long_name = "time interval (days)"

        self.lon,self.lat = {},{}
        self.depth = {}
        for grid in grids:
            mesh = meshes[grid]
            if 'x' in dims and 'y' in dims:
                lon ='glam%s' % grid
                lonNd = f.createVariable(lon,np.float32,htuple)
                Ndname = '%s-longitude' % grid
                lonNd.standard_name  = Ndname
                lonNd.long_name  = Ndname
                self.lon[grid] = lon
                lonNd.units  = 'degrees_east'
                lonNd._CoordinateAxisType = "Lon"
                lonNd[...] = mesh.glam[...]#[1:-1,1:-1]

                if 'glambnds' in mesh._fields:
                    lon_bnds = '%s_bounds' % lon
                    lonNd.bounds =lon_bnds
                    lonbndsNd = f.createVariable(lon_bnds,np.float32,htuple + vtuple)
                    lonbndsNd[...] = mesh.glambnds[...]

                lat = 'gphi%s' % grid
                latNd = f.createVariable(lat,np.float32,('y','x'))
                Ndname = '%s-latitude' % grid
                latNd.standard_name  = Ndname
                latNd.long_name  = Ndname
                self.lat[grid] = lat
                latNd.units  = 'degrees_north'
                latNd._CoordinateAxisType = "Lat"
                latNd[...] = mesh.gphi[...]#[1:-1,1:-1]
                if 'gphibnds' in mesh._fields:
                    lat_bnds = '%s_bounds' % lat
                    latNd.bounds =lat_bnds
                    latbndsNd = f.createVariable(lat_bnds,np.float32,htuple + vtuple)
                    latbndsNd[...] = mesh.gphibnds[...]

                if 'z' in dims:
                    depth = 'gdep%s' % grid
                    depthNd = f.createVariable(depth,np.float32,('z',) + htuple)
                    Ndname = '%s-depth' % grid
                    depthNd.standard_name  = Ndname
                    depthNd.long_name  = Ndname
                    self.depth[grid] = depth
                    depthNd.units  = 'm'
                    depthNd.positive = 'down'
                    depthNd._CoordinateAxisType = "Depth"
                    if mesh.gdep.ndim==1:
                        for j in range(ny):
                            for i in range(nx):
                                depthNd[:,j,i] = mesh.gdep[...]#[1:-1,1:-1]
                    elif mesh.gdep.ndim==3:
                        depthNd[...] = mesh.gdep[...]#[1:-1,1:-1]
                    else:
                        sys.exit('mesh.gdep has shape',mesh.gdep.shape)
            elif 'y' in dims:
                lat = 'gphi%s' % grid
                latNd = f.createVariable(lat,np.float32,('y'))
                Ndname = '%s-latitude' % grid
                latNd.standard_name  = Ndname
                latNd.long_name  = 'northernmost %s' % Ndname
                self.lat[grid] = lat
                latNd.units  = 'degrees_north'
                latNd._CoordinateAxisType = "Lat"
                latNd[...] = mesh.gphi[...].max(-1)#[1:-1,1:-1]
                if 'z' in dims:
                    depth = 'gdep%s' % grid
                    depthNd = f.createVariable(depth,np.float32,('z',) )
                    Ndname = '%s-depth' % grid
                    depthNd.standard_name  = Ndname
                    depthNd.long_name  = Ndname
                    self.depth[grid] = depth
                    depthNd.units  = 'm'
                    depthNd.positive = 'down'
                    depthNd._CoordinateAxisType = "Depth"
                    depthNd[...] = mesh.gdep_0[...]
            #print(f.variables.keys())


    def set_time_index(self,time_index_name):
        if time_index_name == 'years':
            monthlen = np.array([0.,31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.])
            months = ['m%02i'% m for m in np.arange(12)+1] + ['y01']
            endmonths = monthlen.cumsum()
            midmonths = (.5*(endmonths[:-1] + endmonths[1:])).tolist()
            midyear = 365./2.
            midmonths.append(midyear)
            monthlen[:-1] = monthlen[1:]
            monthlen[-1] = 365.
            self.midmonths = dict(zip(months,midmonths))
            self.lengthmonths = dict(zip(months,monthlen))

    def set_tracers(self,tracdict,zlib=False):
        # if fill_value==None:fill_value = netCDF4.default_fillvals['f4']
        self.PNds ={}
        for tracer in tracdict.values():
            trD = tracer.__dict__
            if tracer._FillValue is None:
                fill_value = netCDF4.default_fillvals[tracer.nos.dtype.str[-2:]]
            else:
                fill_value = tracer._FillValue
            #print(tracer.dimensions)
            #print(self.f.variables.keys())
            # list_dims = list(tracer.dimensions)
            # if list_dims[0] != 't':
            #     list_dims[0] = 't'
            #     tracer.dimensions = tuple(list_dims)
            #     print(tracer.dimensions)
            # list(set(tracer.dimensions) & set(self.f.variables.keys()))
            Nd = self.f.createVariable(nemo_mean_names.get(tracer.name,tracer.name),
                                       tracer.nos.dtype,tracer.dimensions,
                                       fill_value=fill_value, zlib=zlib)
            for att in {'coordinates','units','long_name','standard_name',
                        'short_name','volume','area','online_operation'} & set(trD.keys()):
                Nd.__setattr__(att,trD[att])
            tracer.Nd = Nd
            tracer.Nd_lims = {}
            for lim_att in {'median','min','max','var'} & set(trD.keys()):
                Nd = self.f.createVariable(nemo_mean_names.get(tracer.name,tracer.name) +
                                           '_%s' % lim_att,tracer.nos.dtype,tracer.dimensions[:-2],
                                       fill_value=fill_value, zlib=zlib)
                for att in {'units'} & set(trD.keys()):
                    Nd.__setattr__(att,trD[att])
                for att in {'long_name','short_name','standard_name'} & set(trD.keys()):
                    Nd.__setattr__(att,trD[att]+'_%s' % lim_att)
                tracer.Nd_lims[lim_att] = Nd


    def __call__(self,tdict,year=None,month=None,day=None):
        l = self.TC
        self.TCNd[l] = l
        self.TC += 1

        # Unsatisfactory fix
        if year is not None and day is None:
            if month is None:
                month = 'y01'
            else:
                month = 'm%02i'% month

            dtime = self.lengthmonths[month]
            yearfrac = self.midmonths[month]/365.
            time_index_pp5 = np.float32(year + yearfrac)
        else:
            dtime = 5.
            time_index_pp5 = day

        # print 'month/year =',month,dtime
        # print self.lengthmonths
        self.time_indexNd[l] = time_index_pp5
        self.dtime_indexNd[l] = dtime

        #print( tdict.keys())
        for tracer in tdict.values():
            tracer.Nd[l,...] = tracer.nos[...].copy()
            trD = tracer.__dict__
            for att in {'median','min','max','var'} & set(trD.keys()):
                tracer.Nd_lims[att][l,...] = getattr(tracer,att)

    def close(self):
        self.f.close()


class Ddy_sigma(object):
    grid = 'v'
    def __init__(self,name,liked,meshes=None,**kwargs):#,assocd,**kwargs):
        self.grid = self.__class__.grid
        self.mask = ~(meshes[self.grid].mask.astype(bool))
        self.dtype = liked.nos.data.dtype
        self._FillValue = liked._FillValue
        if self._FillValue is None:
            self._FillValue, = np.zeros([1],dtype=liked.nos.dtype)
        self.data = data_like(liked,name)
        self.data.grid = self.grid
        self.data.coordinates = change_str(self.data.coordinates,
                                glamt=u'glam%s' % self.data.grid,
                                gphit=u'gphi%s' % self.data.grid)
        self.describe(**kwargs)

    def working(self,meshes):
        self.Tmask = ~(meshes['t'].mask.astype(bool))

        if self.dtype == np.float32:
            self.sigma_n = eos.sigma_n4
        elif self.dtype == np.float64:
            self.sigma_n = eos.sigma_n8

    def describe(self,**kwargs):
        self.data.long_name = 'meridional sig0 gradient'
        self.data.standard_name = 'd/dy sig0'
        self.data.units = '10^-6 kg m^-3'

    def calc(self,Td,Sd,meshes):
        ddysigma = np.zeros_like(Td.nos.data)
        nz,ny,nx = Td.nos.shape
        for k in range(nz):
            T,S = [x.nos.data[k,...].ravel() for x in (Td, Sd)]
            for j in range(ny-1):
                sig0 = self.sigma_n(self._FillValue,self.Tmask[k,...].ravel(),T,S,0.).reshape([ny,nx])
                ddysigma[k,j,:] = (sig0[j+1,:] - sig0[j,:])*1.e6/meshes['v'].e2[j,:]

        self.data.nos = ma.masked_where(self.mask,ddysigma)

class Ddy_LSPV(Ddy_sigma):
    grid = 'vw'

    def describe(self,**kwargs):
        self.data.long_name = 'meridional LSPV gradient'
        self.data.standard_name = 'd/dy lspv'
        self.data.units = '10^-12 kg m^-5s^-1'

    def calc(self,lspvd,meshes):
        ddy_lspv = np.zeros_like(lspvd.nos.data)
        nz,ny,nx = ddy_lspv.shape
        for k in range(nz):
            for j in range(ny-1):
                ddy_lspv[k,j,:] = (lspvd.nos[k,j+1,:] - lspvd.nos[k,j,:])*1.e5/meshes['vw'].e2[j,:]

        self.data.nos = ma.masked_where(self.mask,ddy_lspv)

def put_z_inner(x):
    nz,ny,nx = x.shape
    return (x.reshape(nz,-1).T.copy()).reshape(ny,nx,nz)


class Rho(object):
    def __init__(self,name,liked, neos=2, **kwargs):#,assocd,**kwargs):
        eos.eos_init(neos)
        self.neos = neos
        self.mask = liked.nos.mask
        self.dtype = liked.nos.data.dtype
        self.shape = self.mask.shape
        self._FillValue = liked._FillValue
        if self._FillValue is None:
            self._FillValue, = np.zeros([1],dtype=liked.nos.dtype)

        self.data = data_like(liked,name)
        self.describe(**kwargs)

    def working(self,meshes):
        if self.dtype == np.float32:
            self.eos_insitu_m = eos.eos_insitu4_m
        elif self.dtype == np.float64:
            self.eos_insitu_m = eos.eos_insitu4_m

        self.depth = meshes['t'].gdep.ravel()

    def describe(self, **kwargs):
        self.data.long_name = 'in-situ density - 1000 using NEMO EOS'
        self.data.standard_name = 'in-situ density'
        self.data.units = 'kg/m^3'

    def calc(self,Td,Sd):
        T,S = [x.nos.data.ravel() for x in (Td, Sd)]
        p = self.depth
        self.data.nos = ma.masked_where(self.mask,self.eos_insitu_m(self._FillValue,self.mask.ravel(),T,S,p).reshape(self.shape))

class BoussinesqR(object):
    """
    Instead of steric anomaly uses
    r_B = rho_in_situ - rho(z,T0,S0) + rho(0,T0,S0)
    appropriate for models with boussinesq equation of state
    """
    def __init__(self,depth,fillvalue,mask,dtype=np.float64,T0=0.,S0=35.,depth_km=0.,neos=2):
        self.dtype = dtype
        self.zt = depth.ravel()
        print( 'in BoussinesqR shape= ',depth.shape)
        self.shape = depth.shape
        eos.eos_init(neos)
        if dtype == np.float32:
            self.eos_insitu_m = eos.eos_insitu4_m
            self.eos_insitu0_m = eos.eos_insitu04_m
        elif self.dtype == np.float64:
            self.eos_insitu_m = eos.eos_insitu4_m
            self.eos_insitu0_m = eos.eos_insitu04_m
        self.T0 = T0
        self.S0 = S0
        self.depth_km = depth_km
        self.fillvalue = fillvalue
        self.mask = mask.ravel()

    def calculate_drho0(self,T0=None,S0=None):
        if T0 is None:
            T0 = self.T0
        else:
            self.T0 = T0
        if S0 is None:
            S0 = self.S0
        else:
            self.S0 = S0

        self.drho0 = self.eos_insitu0_m(self.fillvalue,self.mask,T0,S0,self.depth_km,self.zt).reshape(self.shape)

    def calculate_rho(self,Tz,Sz,zt=None):
        if zt is None:
            ztr = self.zt
        else:
            ztr = zt.ravel()
        self.rho = self.eos_insitu_m(self.fillvalue,self.mask,Tz.ravel(),Sz.ravel(),ztr).reshape(self.shape)

    def __call__(self):
        """
        Outputs r_B in MKS units m^3/kg
        """
        return (self.rho - self.drho0).reshape(self.shape)

class Montgomery(Rho):
    def working(self,meshes,Td,T0=0,S0=35.,depth_km=0.,deltaT=None,deltaS=None):
        Td.nos = ma.masked_where(np.isnan(Td.nos), Td.nos)
        e3w = put_z_inner(meshes['w'].e3)
        e3t = put_z_inner(meshes['t'].e3)
        depth = put_z_inner(meshes['t'].gdep)
        self.Tmask, = [put_z_inner(x) for x in (Td.nos.mask,)]
        self.kmt = (~(self.Tmask)).astype(np.int32).sum(-1).astype(np.int32)
        self.Falsemask = self.kmt < 0

        # interpolate4, interpolate8, mginterpolate4, mginterpolate8,siginterpolate4
        print( 'in Montgomery dtype is',self.dtype)
        if self.dtype == np.float32:
            self.interpolate = interp.interpolate4
            self.mginterpolate = interp.mginterpolate4
            self.siginterpolate = interp.siginterpolate4
        elif self.dtype == np.float64:
            self.interpolate = interp.interpolate8
            self.mginterpolate = interp.mginterpolate8
            self.siginterpolate = interp.siginterpolate8

        self.e3w = e3w
        self.e3t = e3t
        self.depth = depth

        self.rb = BoussinesqR(self.depth,self._FillValue,self.Tmask,self.dtype,T0=T0,S0=S0,depth_km=depth_km)

        self.depth_km = depth_km
        if deltaT is not None:
            self.deltaT = deltaT
        else:
            self.deltaT = 0.5

        if deltaS is not None:
            self.deltaS = deltaS
        else:
            self.deltaS = 0.2

        self.grav = 9.81
        self.rrho0 = 1./1026.
        self.first_time_level = True

    def describe(self,d0=27.):
        self.data.long_name = 'Boussinesq Montgomery function on r_B = %5.2f' % d0
        self.data.standard_name = 'Montgomery function'
        self.data.units = 'm^2s^-2'
        self.d0 = d0

    def calc(self,sshd,Td,Sd,iterate_TS0=False):
        d0 = self.d0
        depth = self.depth
        t0 = time.time()
        T,S = [put_z_inner(x) for x in (Td.nos.data,Sd.nos.data)]
        t1 = time.time()
        print( 'time to reverse T & S is',t1-t0)
        self.rb.calculate_rho(T,S)
        t0 = time.time()
        print( 'time to calculate rho is',t0-t1)

        JM,IM = self.kmt.shape
        T0,S0 = self.rb.T0,self.rb.S0

        sea_indices = self.kmt.nonzero()
        if len(sea_indices[0]) == 0:
            sys.exit('no sea points in domain')
        t0 = time.time()
        while True:
            t2 = time.time()
            if self.first:
                self.rb.calculate_drho0(T0,S0)
            t4 = time.time()
            print( 'time to calculate drho0 is',t4-t2)
            # d3 = self.rb()
            t3 = time.time()
            print( 'time to calculate d is',t3-t4)
            print( self.kmt.min())
            k_below_s,r_above_s,T_s,S_s,outcropmask,groundmask = \
              [x.T for x in self.interpolate(self.kmt.T,T.T,S.T,self.rb.rho.T,self.rb.drho0.T,d0)]

            outcropmask,groundmask = outcropmask.astype(bool),groundmask.astype(bool)
            t2 = time.time()
            print( 'time to calculate k_below is',t2-t3)
            d0mask = self.mask + groundmask + outcropmask
            active = np.logical_not(d0mask)
            if active.max() == False:
                sys.exit('isopycnal always outcrops or grounds in domain')
            T_act =  T_s[active].ravel()
            Tmin_s,Tmax_s,Tbar = T_act.min(),T_act.max(),np.median(T_act)
            S_act =  S_s[active].ravel()
            Smin_s,Smax_s,Sbar = S_act.min(),S_act.max(),np.median(S_act)
            t3 = time.time()
            print( 'time to calculate medians is',t3-t2)

            print( 'on density=',d0,'Tmin Tmax Tmedian=',Tmin_s,Tmax_s,Tbar)
            print( 'on density=',d0,'Smin Smax Smedian=',Smin_s,Smax_s,Sbar)
            print( 'on density=',d0,'T0=',T0, 'S0=', S0)

           # if first call, iterate; then just take already set value unless specify iterate = True
            if (TS0_iterate is 'none' or TS0_iterate is 'first' and not self.first_time_level) or \
               ( abs(T0 - Tbar) < self.deltaT and abs(S0 - Sbar) < self.deltaS ):
                break
            else:
                T0,S0 = Tbar,Sbar
        t1 = time.time()
        print( 'time to find surface is',t1-t0)

        z_s,Mg = [x.T for x in self.mginterpolate(self.kmt.T,T.T,S.T,self.rb.rho.T,self.rb.drho0.T,sshd.nos.data.T,self.e3w.T,
                                                  self.depth.T,k_below_s.T,r_above_s.T,active.T,self.depth_km,d0)]
        t0 = time.time()
        print( 'time to calculate Montgomery',t0-t1)

        self.active = active
        self.data.nos = ma.masked_where(d0mask,Mg)
        self.z_s = ma.masked_where(d0mask,z_s)
        self.incrop_s = groundmask
        self.outcrop_s = outcropmask
        self.z_median_km = ma.median(self.z_s)*1.e-3
        self.sig_s,self.sig_s_zmed = [ma.masked_where(d0mask,x.T) for x in self.siginterpolate(T.T,S.T,
                                           k_below_s.T,r_above_s.T,active.T, self.depth_km,self.z_median_km)]
        self.T_s = ma.masked_where(d0mask,T_s)
        self.T_s_lims = [Tmin_s,Tbar,Tmax_s]
        self.k_below_s = ma.masked_where(self.mask,k_below_s.astype(np.int8))
        self.r_above_s = ma.masked_where(d0mask,r_above_s.astype(np.float32))
        self.S_s = ma.masked_where(d0mask,S_s)
        self.S_s_lims = [Smin_s,Sbar,Smax_s]

        self.first_time_level = False


class Z_s(Rho):
    def describe(self,montg_instance=None):
        self.data.long_name = 'Depth of r_B = %5.2f' % montg_instance.d0
        self.data.standard_name = 'Layer depth'
        self.data.units = 'm'

    def calc(self,montg_instance):
        self.data.nos = montg_instance.z_s
        self.setlims()

    def setlims(self):
        self.data.min, self.data.median, self.data.max = \
          self.data.nos.min(), ma.median(self.data.nos), self.data.nos.max()

class Outcrop_s(object):
    def __init__(self,name,liked,**kwargs):#,assocd,**kwargs):
        self.mask = liked.nos.mask
        self.dtype = np.uint8
        self.shape = self.mask.shape
        self._FillValue = 0# -127

        self.data = data_like(liked,name)
        self.data._FillValue = self._FillValue
        self.describe(**kwargs)
    def describe(self,montg_instance=None):
        self.data.long_name = 'Outcrop region for  r_B = %5.2f' % montg_instance.d0
        self.data.standard_name = 'Outcrop'
        self.data.units = '#'

    def calc(self,montg_instance):
        self.data.nos = montg_instance.outcrop_s.astype(np.uint8)
        #self.setlims()

class Incrop_s(Outcrop_s):

    def describe(self,montg_instance=None):
        self.data.long_name = 'Incrop region for  r_B = %5.2f' % montg_instance.d0
        self.data.standard_name = 'Incrop'
        self.data.units = '#'

    def calc(self,montg_instance):
        self.data.nos = montg_instance.incrop_s.astype(np.uint8)
        #self.setlims()

class K_below_s(object):
    def __init__(self,name,liked,**kwargs):#,assocd,**kwargs):
        self.mask = liked.nos.mask
        self.dtype = np.int8
        self.shape = self.mask.shape
        self._FillValue = -127

        self.data = data_like(liked,name)
        self.data._FillValue = self._FillValue
        self.describe(**kwargs)

    def describe(self,montg_instance=None):
        self.data.long_name = 'Fortran k index just above surface r_B = %5.2f' % montg_instance.d0
        self.data.standard_name = 'Layer index'
        self.data.units = '#'

    def calc(self,montg_instance):
        self.data.nos = montg_instance.k_below_s


class Passive_s(Z_s):
    def describe(self,montg_instance=None):
        if self.data.name=='PWT_s':
            long_trname = 'Pacific Water tracer'
            standard_name = 'Pacific Water tracer'
            self.data.units = '#'

        self.data.long_name = '%s  on %s' % (long_trname, montg_instance.d0)
        self.data.standard_name = '%s  on %s' % (trname, montg_instance.d0)

    def calc(self,tr, montg_instance):
    #  use method from montgomery instance, which has previously output k_lower,r_upper for w-grid as well
        p_s = tracer_interpolate(tr.nos,montg_instance.k_below_s,montg_instance.r_above_s,montg_instance.active)
        self.data.nos = ma.masked_where(~montg_instance.active, p_s)
        self.setlims()


class Sigma0_s(Z_s):
    def describe(self,montg_instance=None):
        self.data.long_name = 'Potential density on r_B = %5.2f' % montg_instance.d0
        self.data.standard_name = 'Layer sigma0'
        self.data.units = 'kg/m^3'

    def calc(self,montg_instance):
        self.data.nos = montg_instance.sig_s
        self.setlims()

    def setlims(self):
        self.data.min, self.data.var, self.data.max = \
          self.data.nos.min(), ma.var(self.data.nos), self.data.nos.max()


class SigmaMedian_s(Sigma0_s):
    def describe(self,montg_instance=None):
        self.data.standard_name = 'Layer sigma'
        self.data.units = 'kg/m^3'

    def calc(self,montg_instance):
        self.data.nos = montg_instance.sig_s_zmed
        self.data.long_name = \
          'Potential density ref to %5.2f km on r_B = %5.2f' %(montg_instance.z_median_km, montg_instance.d0)
        self.setlims()


class T_s(Rho):
    def describe(self,montg_instance=None):
        self.data.long_name = 'Temperature on r_B = %5.2f' % montg_instance.d0
        self.data.standard_name = 'Layer temperature'
        self.data.units = 'deg C'

    def calc(self,montg_instance):
        self.data.nos = montg_instance.T_s
        lims = montg_instance.T_s_lims
        self.setlims(lims)

    def setlims(self,lims):
        self.data.min, self.data.median, self.data.max = lims

class S_s(T_s):
    def describe(self,montg_instance=None):
        self.data.long_name = 'Salinity on r_B = %5.2f' % montg_instance.d0
        self.data.standard_name = 'Layer salinity'
        self.data.units = 'psu'

    def calc(self,montg_instance):
        self.data.nos = montg_instance.S_s
        lims = montg_instance.S_s_lims
        self.setlims(lims)



class BottomPressure(Rho):
    def working(self,meshes,Td):
        if self.dtype == np.float32:
            self.bp_mn = bp.bp.bp_mn4
        elif self.dtype == np.float64:
            self.bp_mn = bp.bp.bp_mn8

        self.depth = meshes['t'].gdep.ravel()
        e3w = meshes['w'].e3.T.copy()
        kmt = (~(Td.nos.mask)).astype(np.int32).sum(0).astype(np.int32)
        # print e3w.dtype,'\n',e3w.flags,'\n',kmt.dtype,'\n',kmt.flags
        bp.bp.setup_bp_mn(self._FillValue,kmt.T,e3w.T)

    def describe(self):
        self.data.long_name = 'bottom pressure anomaly from rh=1026/ssh=0 using NEMO EOS'
        self.data.standard_name = 'bottom pressure'
        self.data.units = 'kPa'

    def calc(self,sshd,rhod=None,Td=None,Sd=None):
        rho, = [x.T.copy() for x in (rhod.nos.data,)]
        # print rho.dtype,'\n',rho.flags,'\n',sshd.nos.data.dtype,'\n',sshd.nos.data.flags
        self.data.nos = ma.masked_where(self.mask,self.bp_mn(rho.T,sshd.nos.data.T).T)


class TotalHeat(Rho):
    def working(self):
        self.cp = 4.e3

    def describe(self):
        self.data.long_name = 'total downward heat H_in - (E-P)*c_p*SST'
        self.data.standard_name = 'total heat in'
        self.data.units = 'W/m^2'

    def calc(self,Hind,EmPd,sstd):
        self.data.nos = Hind.nos - self.cp*EmPd.nos*sstd.nos

class WaterIce(Rho):
    def working(self):
        self.rh_ice = 900.
        self.rrho0 = 1./1026.

    def describe(self):
        self.data.long_name = 'Melted ice height aice*hice*rh_ice/rho0'
        self.data.standard_name = 'Melted ice'
        self.data.units = 'm'

    def calc(self,aiced,hiced):
        self.data.nos = self.rrho0*(1.-aiced.nos)*self.rh_ice*hiced.nos

class WaterSnow(Rho):
    def working(self):
        self.rh_snow = 330.
        self.rrho0 = 1./1026.

    def describe(self):
        self.data.long_name = 'Melted snow height aice*hsnow*rh_snow/rho0'
        self.data.standard_name = 'Melted snow'
        self.data.units = 'm'

    def calc(self,aiced,hsnowd):
        self.data.nos = self.rrho0*(1.-aiced.nos)*self.rh_snow*hsnowd.nos


class ExplicitHeat(TotalHeat):
    def describe(self):
        self.data.long_name = 'explicit downward heat H_in'
        self.data.standard_name = 'explicit heat in'
        self.data.units = 'W/m^2'

    def calc(self,Hind):
        self.data.nos = Hind.nos

class SaltIn(Rho):
    def describe(self):
        self._FillValue=None
        self.data.long_name = 'salt from ice melting/freezing into ocean:[((E-P)_s - (E-P)]*(SS0/rho0)'
        self.data.standard_name = 'salt in'
        self.data.units = 'kg/(m^2 s)'

    def calc(self,EmPd,EmPsd,sssd):
        # self.data.nos = (EmPsd.nos - EmPd.nos)*(35./1026.)
        self.data.nos = (EmPsd.nos - EmPd.nos)*(sssd.nos/1026.)

class WaterIn(Rho):
    def describe(self):
        self.data.long_name = 'total water input into ocean'
        self.data.standard_name = 'water input'
        self.data.units = 'kg/(m^2 s)'

    def calc(self,EmPd):
        self.data.nos = -EmPd.nos


class BuoyancyIn(Rho):
    def describe(self):
        self.data.long_name = 'total buoyancy (Boussinesq mass) input into ocean'
        self.data.standard_name = 'mass input'
        self.data.units = '10^-6 kg/(m^2 s)'

    def working(self):
        if self.dtype == np.float32:
            self.eos_rab_ref_m = eos.eos_rab_ref4_m
        # elif self.dtype == np.float64:
        #     self.rho_mn = eos.eos_insitu4_m

    def calc(self,EmPsd,Hind,sstd,sssd):
        rcp = 4.e3
        salt =  EmPsd.nos*(sssd.nos/1026.)
        sst,sss = [x.nos.data.ravel() for x in (sstd, sssd)]
        alpha, beta = self.eos_rab_ref_m(self._FillValue,self.mask.ravel(),sst,sss,0.)
        self.data.nos = 1.e6*(beta.reshape(self.shape)*salt*1026. -
                         alpha.reshape(self.shape)*Hind.nos/rcp)


class Glob3Av(object):
    avtype = 'Global mean'

    def __init__(self,name,Td):
        self.data = Data(name)
        TdD = Td.__dict__
        self.data._FillValue = None
        for att in {'grid','units','online_operation'} & set(TdD.keys()):
            self.data.__setattr__(att,TdD[att])
        for att in {'long_name','standard_name','short_name'} & set(TdD.keys()):
            self.data.__setattr__(att,'%s  %s' % (self.__class__.avtype,TdD[att]))
        self.get_dimensions()
        print( 'In init  of ',name,' dimensions are ',self.data.dimensions)
        self.describe()

    def get_dimensions(self):
        self.data.dimensions = ('time_counter',)

    def describe(self):
        pass

    def working(self,meshes):
        mesh = meshes[self.data.grid]
        self.dA = mesh.e1*mesh.e2*mesh.maskutil
        self.dV = mesh.mask*mesh.e3*self.dA[None,:,:]
        self.V0 = self.dV.sum()
        self.data.volume = self.V0

    def calc(self,Td,sshd):
        dV = self.dV.copy()
        dV0 = sshd.nos*self.dA
        dV[0,...] += dV0
        V = self.V0 + dV0.sum()
        self.data.nos = np.dot(dV.ravel(),Td.nos.data.ravel())/V


class Glob2Av(Glob3Av):
    def working(self,meshes):
        mesh = meshes[self.data.grid]
        self.dA = mesh.e1*mesh.e2*mesh.maskutil
        self.A = self.dA.sum()
        self.data.area = self.A

    def calc(self,Td):
        self.data.nos = np.dot(self.dA.ravel(),Td.nos.data.ravel())/self.A


class Glob3Sum(Glob3Av):
    avtype = 'Global 3d Sum'

    def calc(self,Td,sshd):
        dV = self.dV.copy()
        dV0 = sshd.nos*self.dA
        dV[0,...] += dV0
        self.data.nos = np.dot(dV.ravel(),Td.nos.data.ravel())

class Glob2Sum(Glob2Av):
    avtype = 'Global 2d Sum'

    def calc(self,Td):
        self.data.nos = np.dot(self.dA.ravel(),Td.nos.data.ravel())

class Glob3Heat(Glob3Sum):
    def describe(self):
        self.data.long_name = 'ocean heat content within domain'
        self.data.standard_name = 'heat content'
        self.data.short_name = 'heat content'
        self.data.units = 'J'

    def working(self):
        cp = 4.e3
        rho0 = 1026.
        self.cprho = cp*rho0

    def calc(self,Tsum):
        self.data.nos = Tsum.nos*self.cprho


class Glob3FW(Glob3Sum):
    def describe(self):
        self.data.units = 'm^3'

    def working(self,S0=35.):
        self.S0 = S0
        self.data.long_name = 'ocean FW content within domain, ref S= %2.1f' % S0
        self.data.standard_name = 'FW content'
        self.data.short_name = 'FW content'

    def calc(self,Ssum,sshsum):
        V = Ssum.volume + sshsum.nos
        self.data.nos = V - Ssum.nos/self.S0


class Zon3Av(object):
    def __init__(self,name,liked, **kwargs):#,assocd,**kwargs):
        self.mask = liked.nos.mask
        self.dtype = liked.nos.data.dtype
        self.shape = self.mask.shape
        self._FillValue = liked._FillValue
        if self._FillValue is None:
            self._FillValue, = np.zeros([1],dtype=liked.nos.dtype)

        self.data = data_like(liked,name)
        dimensions_list = list(self.data.dimensions)
        dimensions_list.remove('x')
        self.data.dimensions = tuple(dimensions_list)
        coordinates_list = self.data.coordinates.split()
        coordinates_list = [c for c in coordinates_list if 'glam' not in c]
        depth_name = 'gdep%s' % self.data.grid
        for i,c in enumerate(coordinates_list):
            if 'dep' in c:
                coordinates_list[i] = depth_name
                break
        else:
            coordinates_list.append(depth_name)
        self.data.coordinates = ' '.join(coordinates_list)
        trD = liked.__dict__
        self.data.units = liked.units
        self.data.standard_name = 'zonal av %s' % trD.get('standard_name', self.data.name)
        self.data.long_name = 'zonally averaged %s' % trD.get('long_name', self.data.standard_name)
        self.describe(**kwargs)
        # print self.data.__dict__


    def describe(self,**kwargs):
        pass

    def working(self,meshes):
        mesh = meshes[self.data.grid]
        # self.dA = mesh.e1*mesh.e2*mesh.maskutil
        # self.dV = mesh.mask*mesh.e3*self.dA[None,:,:]
        # self.dA = ma.masked_where(self.mask[0,...],mesh.e1*mesh.e2*mesh.maskutil).filled(0.)
        # self.dV = ma.masked_where(self.mask,mesh.mask*mesh.e3*self.dA[None,:,:]).filled(0.)
        self.dA = mesh.e1*mesh.e2*mesh.maskutil
        dV = mesh.mask*mesh.e3*self.dA[None,:,:]
        dV0 = dV[0,...]
        V = dV.sum(-1)
        self.dV0 = dV[0,...]
        self.mask = V<=1.e-5
        V[self.mask] = 1.e-5
        self.V0 = V[0,...]
        self.weight = dV/V[:,:,None]
        # print self.weight[10,259,:].sum(-1)
        print( 'calculated mask')

    def calc(self,Td,sshd):
        weight = self.weight
        nz,ny,nx = weight.shape
        zonav = np.empty([nz,ny],dtype=Td.nos.data.dtype)

        if sshd is None or self.data.grid not in ['w','t']:
            for k in range(nz):
                for j in range(ny):
                    zonav[k,j] = np.dot(weight[k,j,:],Td.nos.data[k,j,:])
            k,j = 10,30
            # print 'zonav[10,30] =', zonav[k,j]
        else:
            dV0 = sshd.nos*self.dA
            V0 = dV0.sum(-1) + self.V0[...]
            weight0 = (dV0 + self.dV0)/V0[:,None]
            k = 0
            for j in range(ny):
                zonav[k,j] = np.dot(weight0[j,:],Td.nos.data[k,j,:])

            for k in range(2,nz):
                for j in range(ny):
                    zonav[k,j] = np.dot(weight[k,j,:],Td.nos.data[k,j,:])

        self.data.nos = ma.masked_where(self.mask,zonav)
        # print 'dimensions',self.data.dimensions

@jit
def jacobian2(quarter_rarea,r1,r2):
    nz,ny,nx = r1.shape
    jacobian = np.zeros([nz,ny,nx],dtype=r1.dtype)
    jacobian_angle = np.zeros([nz,ny,nx],dtype=r1.dtype)
    dyr1 = np.zeros([nx],dtype=jacobian.dtype)
    dyr2 = np.zeros([nx],dtype=jacobian.dtype)
    dxr1 = np.zeros([2,nx-1],dtype=jacobian.dtype)
    dxr2 = np.zeros([2,nx-1],dtype=jacobian.dtype)

    rad2deg = 180./math.pi

    for k in range(nz):
        r1k0,r2k0 = r1[k,0,:],r2[k,0,:]
        dxr1[0,:] = r1k0[1:] - r1k0[:-1]
        dxr2[0,:] = r2k0[1:] - r2k0[:-1]
        l = 1
        for j in range(ny-1):
            l0 = 1 - l
            r1k1,r2k1 = r1[k,j+1,:],r2[k,j+1,:]
            dxr1[l,:] = r1k1[1:] - r1k1[:-1]
            dxr2[l,:] = r2k1[1:] - r2k1[:-1]

            dyr1[:] = r1k1 - r1k0
            dyr2[:] = r2k1 - r2k0

            for i in range(nx-1):
                jacobian[k,j,i] = quarter_rarea[j,i]*(
                                       dxr1[l0,i]*dyr2[i] - dxr2[l0,i]*dyr1[i]
                                     + dxr1[l0,i]*dyr2[i+1] - dxr2[l0,i]*dyr1[i+1]
                                     + dxr1[l,i] *dyr2[i]   - dxr2[l,i]*dyr1[i]
                                     + dxr1[l,i] *dyr2[i+1] - dxr2[l,i]*dyr1[i+1]
                                                     )
                gradmag = 2.*quarter_rarea[j,i]*math.sqrt(
                                        ( dyr1[i]**2  +  dyr1[i+1]**2
                                      +  dxr1[l,i]**2 + dxr1[l0,i]**2 )
                                        * ( dyr2[i]**2  +  dyr2[i+1]**2
                                      +  dxr2[l,i]**2 + dxr2[l0,i]**2 )
                                                         )
                jacobian_angle[k,j,i] = rad2deg*math.asin(
                    max(min(jacobian[k,j,i]/(gradmag + 1.e-30),1.),-1.)
                                                         )
            r1k0,r2k0 = r1k1,r2k1
            l = l0
    return jacobian, jacobian_angle


class JacobianDepth(Ddy_sigma):
    grid = 'f'

    def describe(self,**kwargs):
        self.data.long_name = 'Jacobian of tracers at constant depth'
        self.data.standard_name = 'Jacobian_depth'
        self.data.units = 'tracers /m^3'

    def working(self,meshes):
        self.quarter_rarea = .25/(meshes[self.grid].e1*meshes[self.grid].e2)

    def calc(self,r1,r2,meshes):
        t0 = time.time()
        jacobian,self.jacobian_angle = jacobian2(self.quarter_rarea,r1.nos,r2.nos)
        t1 = time.time()
        print( 'calculating jacobian took %f seconds' % (t1-t0))
        # jacobian = jacobian2(self.quarter_rarea,r1.nos,r2.nos)
        # t2 = time.time()
        # print 'calculating jacobian 2nd time took %f seconds' % (t2-t1)
        self.data.nos = ma.masked_where(self.mask,jacobian)

    def old_calc(self,r1,r2,meshes):
        nz,ny,nx = r1.nos.shape
        jacobian = np.zeros([nz,ny,nx],dtype=r1.nos.dtype)
        dyr1 = np.zeros([nx],dtype=jacobian.dtype)
        dyr2 = np.zeros([nx],dtype=jacobian.dtype)
        dxr1 = np.zeros([2,nx-1],dtype=jacobian.dtype)
        dxr2 = np.zeros([2,nx-1],dtype=jacobian.dtype)

        for k in range(nz):
            r1k0,r2k0 = r1.nos[k,0,:],r2.nos[k,0,:]
            dxr1[0,:] = r1k0[1:] - r1k0[:-1]
            dxr2[0,:] = r2k0[1:] - r2k0[:-1]
            l = 1
            for j in range(ny-1):
                l0 = 1 - l
                r1k1,r2k1 = r1.nos[k,j+1,:],r2.nos[k,j+1,:]
                dxr1[l,:] = r1k1[1:] - r1k1[:-1]
                dxr2[l,:] = r2k1[1:] - r2k1[:-1]

                dyr1[:] = r1k1 - r1k0
                dyr2[:] = r2k1 - r2k0

                for i in range(nx-1):
                    jacobian[k,j,i] = self.quarter_rarea[j,i]*(   dxr1[l0,i]*dyr2[i]   - dxr2[l0,i]*dyr1[i]
                                            + dxr1[l0,i]*dyr2[i+1] - dxr2[l0,i]*dyr1[i+1]
                                            + dxr1[l,i] *dyr2[i]   - dxr2[l,i]*dyr1[i]
                                            + dxr1[l,i] *dyr2[i+1] - dxr2[l,i]*dyr1[i+1])
                r1k0,r2k0 = r1k1,r2k1
                l = l0

        self.data.nos = ma.masked_where(self.mask,jacobian)

class JacobianAngle(JacobianDepth):
    def describe(self,**kwargs):
        self.data.long_name = 'Jacobian angle of tracers at constant depth'
        self.data.standard_name = 'Jacobian_angle'
        self.data.units = 'degrees'

    def calc(self,jacobian_instance):
        self.data.nos = ma.masked_where(self.mask,jacobian_instance.jacobian_angle)


if __name__ == '__main__':
    t00 = time.time()

    parser = ArgumentParser(description='Output various NEMO diagnostics')

    parser.add_argument('--meshdir',dest='meshdir',
                        help='name of mesh directory; can be set from environment variable MESHDIR',
                        default=None)
    parser.add_argument('--meshfile',dest='meshfile',
                        help='name of meshfile inside mesh directory')

    parser.add_argument('--infile',dest='infile',
                        help='list of path of  data files to read; can include wildcard expressions',
                        nargs= '*',default=[])
    parser.add_argument('--hlimits',dest='hlimits',
                        help='horizontal limits; required order is ylo, yhi,xlo,xhi',
                        nargs=4,type=int,default=[1,-1,1,-1])

    parser.add_argument('--xlimits',dest='xlimits',
                        help='horizontal limits; required order is xlo,xhi',
                        nargs=2,type=int,default=[1,-1])
    parser.add_argument('--ylimits',dest='ylimits',
                        help='horizontal limits; required order is ylo,yhi',
                        nargs=2,type=int,default=[1,-1])

    parser.add_argument('-t','--tracers', dest='mtracers',help='names of mean tracers to read',
                        nargs= '*',default=[])
    parser.add_argument('--passive_s',dest='passive_s',
                        help='names of output passive tracers on surfaces',
                        nargs='*',default=[])
    parser.add_argument('-x','--xtracers',dest='xtracers',
                        help='names of calculated tracers to output',
                        nargs= '*',default=[])

    parser.add_argument('--density',dest='density',
                        help='layer density for layer output',type=float,default=None)
    parser.add_argument('--depth_km',dest='depth_km',
                        help='reference depth in km for sigma_s',type=float,default=0.)
    parser.add_argument('--TS0',dest='TS0',
                        help='initial guess for T0 and S0 on density layer',
                        nargs=2,type=float,default=None)
    parser.add_argument('--TS0_iterate',dest='TS0_iterate',
                        help='How to iterate T0 and S0 in definition of r_b:\n'
                        'none => use T0 and S0 directly as specified from args.TS0'
                        'first => iterate T0 and S0 when doing first time-level; use these for later levels'
                        'all => iterate T0 and S0 independedently syatying from args.Ts0 on each time-level',
                        choices= ['none','first','all'], default='none')

    parser.add_argument('--neos',dest='neos',type=int,default=2,
                        help='choose EOS: -1=> old Jackett McDougall, 0=> poly EOS-80, 2=> TEOS-10')
    parser.add_argument('--nthreads',dest='nthreads',
                        help='number of threads for EOS; 16 is good for workstations',
                        type=int,default=4)
    parser.add_argument('--dims',dest='dims',help='dimensions in output cdf file',
                        nargs='+',default = ['t','z','y','x'])

    parser.add_argument('-r',dest='runname',help='Run name',default='ORCA1-N403')
    parser.add_argument('--xroot',dest='xroots',help='Extra root directories',nargs= '*',default=None)
    parser.add_argument('--refresh',dest='refresh',help='refresh cache?',
                        action='store_true',default=False)
    parser.add_argument('-y',dest='years',help='last and first years',
                        type=int,nargs='+',default=None)
    parser.add_argument('--no_bounds',dest='no_bounds',
                        help='do not output meshbounds',
                        action='store_true',default=False)
    parser.add_argument('--checkmask',dest='checkmask',
                        help='always use mask from mask.nc',
                        action='store_true',default=False)
    parser.add_argument('-y0',dest='year00',help='first year of dataset',
                        type=int,default=1948)
    parser.add_argument('--pass',dest='passno',help='pass number',
                        type=int, nargs='*',default=None)
    parser.add_argument('-m','--months',dest='months',
                        help='months to save',
                        nargs= '*',type=int,default=None)
    parser.add_argument('-d','--days',dest='days',
                        help='days to save',
                        nargs= '*',type=int,default=None)
    parser.add_argument('-o',dest='outdir',help='directory of output file',default='.')
    parser.add_argument('--restarts',dest='rtracers',
                        help='name of restart tracers to read',
                        nargs= '*',default=[])
    args = parser.parse_args()

    eos.set_eos_threads(args.nthreads)
    eos.eos_init(args.neos)
    interp.set_eos_threads(args.nthreads)
    interp.eos_init(args.neos)

    if args.meshdir is None:
        try:
            args.meshdir = os.environ["MESHDIR"]
        except:
            print('meshdir not found from args or environment variable')
            pass

    if not args.infile and args.meshdir is None:
        model_run = findnemo.MultiDir(args.runname,xroots=args.xroots,refresh=args.refresh)

    if args.years is  None:
        years = [None]
        pass_prefix = ''
        passno0 = None
        year0 = None
    else:
        if args.passno is not None:
            passno0 = args.passno[0]
            passno1 = args.passno[-1]
            passes = ''.join([ '%i' % p for p in args.passno])

            year1 = (passno1 - passno0)*60 + args.years[-1]

            pass_prefix = 'P%s' % passes
        else:
            pass_prefix = ''
            passno0 = None
            year1 = args.years[-1]

        year0 = args.years[0]
        nyears = year1 - year0 + 1
        years = np.arange(nyears) + year0

    allmonths = np.arange(12) + 1
    if args.months is None:
        months = [None]
        fs = 'y'
    elif args.months == []:
        months = allmonths
        fs = 'm_all'
    else:
        months = args.months
        fs = 'm' + ''.join([ '%02i' % i for i in args.months])

    if args.days is None:
        days = [None]
        fs = 'y'
    else:
        days = args.days
        fs = 'd' + ''.join([ '%04i' % i for i in args.days])

    filend = '_'
    if args.rtracers:
        filend = '_%s_r' % ''.join(args.rtracers)
    if args.mtracers:
        filend = '_%s_m' % ''.join(args.mtracers)

    if year0 is None:
        fileyrs = ''
    elif year0 == args.years[-1]:
        fileyrs = '%sy%4i%s%s' % (pass_prefix,year0,fs,filend)
    else:
        fileyrs = '%sy%4ito%4i%s%s' % (pass_prefix,year0,args.years[-1],fs,filend)

    tracstr = '_'.join([t.replace('_s','') for t in args.xtracers])
    if args.density is not None:
        tracstr = '%s_%g' % (tracstr,args.density)

    outname = '%s_%s_%s.nc' % (args.runname,fileyrs,tracstr)
    # sys.exit(outname)

    def make_2_slices(hbounds, wideslice=True):

        if wideslice:
            widehbounds = hbounds[:]
            for n,l in enumerate(widehbounds):
                if n % 2 == 0: # lower limits
                    widehbounds[n] -= 1
                else:          # upper limits
                    if l==-1:
                        widehbounds[n] = None
                    else:
                         widehbounds[n] += 1
            hboundslist = [hbounds,widehbounds]
        else:
            hboundslist = [hbounds]

        def make_slice(bounds2d):
            return (slice(bounds2d[0],bounds2d[1],None),slice(bounds2d[2],bounds2d[3],None))

        return [make_slice(hb) for hb in hboundslist]

    hslice,wide_hslice = make_2_slices(args.hlimits)
    t01 = time.time()
    print (f'took {t01-t00} to set args')
    print (f'time at after setting args is {t01-t00}')


    fexttracm = FextTrac('mean')
    fexttracr = FextTrac('restart')
    gridtrac = GridTrac()
    DCDF4.set_gridtrac(gridtrac)
    DCDF4.set_default_slice((0,Ellipsis) + hslice)

    gridtracerd = gridtrac.get_tracdict(args.mtracers + args.rtracers)
    grids = list(gridtracerd.keys())
    grids.sort()

    fexttracerd = fexttracm.get_tracdict(args.mtracers)
    fexttracerr = fexttracr.get_tracdict(args.rtracers)
    fexttracerd.update(fexttracerr)
    fexts = list(fexttracerd.keys())
    fexts.sort()

    inargs = InArgs(args.xtracers,vdict= {'bp':'rho','jacobian_angle':'jacobian_depth','heatsum':'gsT','FWsum':['gsS','gsssh']},
                    start=['zm','gs'],end='bar')
    xgrid = []
    if inargs('bp') or inargs('mont'):
        xgrid.append('w')
    if inargs('ddy_lspv'):
        xgrid.append('vw')
    if inargs('ddysigma'):
        xgrid.append('v')
    if inargs('jacobian_depth'):
        xgrid.append('f')
    print( 'grids are',xgrid+grids)
    gridgen = RealGridStuff(grids + xgrid)
    meshkeys = ['glam','gphi','gdep','gdep_0','e3','e1','e2','maskutil','mask']
    if args.meshdir is None:
        meshdir = model_run
    else:
        meshdir = args.meshdir
    meshes = gridgen.find_meshes(meshkeys,meshdir,hslice,
                                 wide_hslice=wide_hslice,meshfile=args.meshfile, meshbnds=not args.no_bounds)
    # print meshes

    # data = {}
    t0 = time.time()
    print (f'took {t0-t01} to set meshes')
    print (f'time at after setting meshes is {t0-t00}')
    tdict = {}
    xdict = {}

    if args.infile:
        args.infile.sort()
        use_infile = True
        print (f'will process files')
        print ('\n'.join(args.infile),'\n')
        infile = args.infile.pop(0)
    else:
        use_infile = False

    for fext in fexts:
        if use_infile:
            pathname = infile[:-3]
            while pathname[-1].isalpha():
                pathname = pathname[:-1]

            pathname += fext + '.nc'
        else:
            pathname = model_run(passno=passno0,years=year0,month=months[0],fext=fext,day5=days[0])


        cdf_file = DCDF4(pathname, checkmask=args.checkmask)
        print(fexttracerd[fext])
        P = cdf_file(fexttracerd[fext], meshes=meshes)
        tdict.update(P)


    t01 = time.time()
    print (f'took {t01-t0} to set output file')
    print (f'time at after setting up output file is {t01-t00}')

    other_dims_dict = {}
    for dataset in tdict.values():
        # print dataset.__dict__
        y = dataset.otherdims
        # print y
        other_dims_dict.update(y)
    other_dims = set(other_dims_dict.keys())

    idict = tdict.copy()

    t0 = time.time()
    print (f'took {t0-t01} to set up dicts')
    print (f'time at after setting up dicts is {t0-t00}')

    if inargs('rho'):
        rho = Rho('rho',tdict['T'], neos=args.neos)
        rho.working(meshes)
        rho.calc(tdict['T'],tdict['S'])
        idict['rho'] = rho.data

    if inargs('TS'):
        TS = BinTS('TS',tdict['T'])
        TS.working(meshes)
        TS.calc(tdict['T'],tdict['S'])
        idict['TS'] = TS.data

    if inargs('ddy_lspv'):
        ddy_lspv = Ddy_LSPV('ddylspv',tdict['lspv'],meshes)
        ddy_lspv.calc(tdict['lspv'],meshes)
        idict['ddy_lspv'] = ddy_lspv.data

    if inargs('ddysigma'):
        ddysigma = Ddy_sigma('ddysigma',tdict['T'],meshes)
        ddysigma.working(meshes)
        ddysigma.calc(tdict['T'],tdict['S'],meshes)
        idict['ddysigma'] = ddysigma.data

    if inargs('jacobian_depth'):
        jacobian_depth = JacobianDepth('jacobian_depth',tdict['T'],meshes)
        jacobian_depth.working(meshes)
        jacobian_depth.calc(tdict['T'],tdict['S'],meshes)
        idict['jacobian_depth'] = jacobian_depth.data

    if inargs('jacobian_angle'):
        jacobian_angle = JacobianAngle('jacobian_angle',tdict['T'],meshes)
        jacobian_angle.calc(jacobian_depth)
        idict['jacobian_angle'] = jacobian_angle.data

    if inargs('bp'):
        bpd = BottomPressure('bp',tdict['ssh'])
        bpd.working(meshes,tdict['T'])
        bpd.calc(tdict['ssh'],rho.data)
        idict['bp'] = bpd.data

    if inargs('mont'):
        mgd = Montgomery('mont',tdict['ssh'],neos=args.neos, d0=args.density)
        T0,S0 = args.TS0
        mgd.working(meshes,tdict['T'],T0=T0,S0=S0,depth_km=args.depth_km)
        mgd.calc(tdict['ssh'],tdict['T'],tdict['S'])
        idict['mont'] = mgd.data

        layerdict = {'k_s':K_below_s,'z_s':Z_s,'sigma_s':Sigma0_s,
                     'sigma_med_s':SigmaMedian_s,'T_s':T_s,'S_s':S_s,
                     'outcrop_s':Outcrop_s, 'incrop_s':Incrop_s}
        instancedict = {}
        for x,y in layerdict.items():
            if inargs(x):
                instance = y(x,tdict['ssh'],montg_instance=mgd)
                instance.calc(mgd)
                instancedict[x] = instance
                idict[x] = instance.data

    passive_s_dict={}
    for x in args.passive_s:
        trname = x[:-2]
        pp = Passive_s(x,tdict['ssh'],montg_instance=mgd)
        pp.calc(tdict[trname], mgd)
        passive_s_dict[x] = pp
        idict[x] = pp.data

    if inargs('Sin'):
        Sind = SaltIn('Sin',tdict['Hin'])
        Sind.calc(tdict['EmP'],tdict['EmPs'],tdict['sss'])
        idict['Sin'] = Sind.data

    if inargs('Bin'):
        Bind = BuoyancyIn('Bin',tdict['Hin'])
        Bind.working()
        Bind.calc(tdict['EmPs'],tdict['Hin'],tdict['sst'],tdict['sss'])
        idict['Bin'] = Bind.data

    if inargs('Win'):
        Wind = WaterIn('Win',tdict['Hin'])
        Wind.calc(tdict['EmP'])
        idict['Win'] = Wind.data

    if inargs('THin'):
        THind = TotalHeat('THin',tdict['Hin'])
        THind.working()
        THind.calc(tdict['Hin'],tdict['EmP'],tdict['sst'])
        idict['THin'] = THind.data

    if inargs('IceW'):
        IceWd = WaterIce('IceW',tdict['aice'])
        IceWd.working()
        IceWd.calc(tdict['aice'],tdict['hice'])
        idict['IceW'] = IceWd.data

    if inargs('SnowW'):
        SnowWd = WaterSnow('SnowW',tdict['aice'])
        SnowWd.working()
        SnowWd.calc(tdict['aice'],tdict['hsnow'])
        idict['SnowW'] = SnowWd.data

    t01 = time.time()
    print (f'took {t01-t0} to calculate variables first')
    print (f'time  after first calculations is {t01-t00}')

    bar2d,bar3d = {},{}

    for xbar in inargs.arguments:
        if xbar.endswith('bar'):
            x = xbar[:-3]
            if nemo_dimensions[x] == 2:
                bar2d[x] = Glob2Av(xbar,idict[x])
                bar2d[x].working(meshes)
                bar2d[x].calc(idict[x])
                idict[xbar] = bar2d[x].data

            elif nemo_dimensions[x] == 3:
                bar3d[x] = Glob3Av(xbar,idict[x])
                bar3d[x].working(meshes)
                bar3d[x].calc(idict[x],tdict['ssh'])
                idict[xbar] = bar3d[x].data

    zm3d = {}
    for zm in inargs.arguments:
        if zm.startswith('zm'):
            x = zm[2:]
            zm3d[x] = Zon3Av(zm,idict[x])
            zm3d[x].working(meshes)
            zm3d[x].calc(idict[x],tdict.get('ssh'))
            idict[zm] = zm3d[x].data

    gsdict = {}
    for gs in inargs.arguments:
        if gs.startswith('gs'):
            x = gs[2:]
            if len(idict[x].nos.shape)==2:
                gsdict[x] = Glob2Sum(gs,idict[x])
            elif len(idict[x].nos.shape)==3:
                gsdict[x] = Glob3Sum(gs,idict[x])
            gsdict[x].working(meshes)
            if len(idict[x].nos.shape)==2:
                gsdict[x].calc(idict[x])
            elif len(idict[x].nos.shape)==3:
                gsdict[x].calc(idict[x],tdict['ssh'])
            idict[gs] = gsdict[x].data


    if inargs('heatsum'):
        G3Heatd = Glob3Heat('heatsum',idict['gsT'])
        G3Heatd.working()
        G3Heatd.calc(idict['gsT'])
        idict['heatsum'] = G3Heatd.data

    if inargs('FWsum'):
        G3FWd = Glob3FW('FWsum',idict['gsS'])
        G3FWd.working()
        G3FWd.calc(idict['gsS'],idict['gsssh'])
        idict['FWsum'] = G3FWd.data

    t0 = time.time()
    print (f'took {t0-t01} to calculate sums, means etc')
    print (f'time  after calculating sums, means etc is {t01-t00}')

    for tname in args.xtracers + args.passive_s:
        # print tname
        if tname in idict.keys():
            # print tname
            xdict[tname] = idict[tname]

    other_dims_dict = {}
    for dataset in xdict.values():
        if hasattr(dataset,'otherdims'):
            y = dataset.otherdims
            other_dims_dict.update(y)
    other_dims = set(other_dims_dict.keys())
    other_dims_dict = {d:other_dims_dict[d] for d in other_dims_dict}


    outgrids = list({xdict[tname].grid for tname in xdict.keys()})

    t01 = time.time()
    print (f'took {t01-t0} to set up output variables')
    print (f'time  after setting up output variables is {t01-t00}')

    threedcdf = Create3DCDF(outgrids,meshes,
                          outpath=pjoin(args.outdir,outname),
                          time_index_name='years',dims=args.dims,
                          other_dims_dict = other_dims_dict,density = args.density)

    threedcdf.set_tracers(xdict,zlib=True)

    t01 = time.time()
    print (f'took {t01-t0} to set up output details into output file')
    print (f'time  after setting output details into output file is {t01-t00}')

    def do_on_file():
            if inargs('rho'):
                rho.calc(tdict['T'],tdict['S'])

            if inargs('TS'):
                TS.calc(tdict['T'],tdict['S'])

            if inargs('ddy_lspv'):
                ddy_lspv.calc(tdict['lspv'],meshes)

            if inargs('ddysigma'):
                ddysigma.calc(tdict['T'],tdict['S'],meshes)

            if inargs('jacobian_depth'):
                jacobian_depth.calc(tdict['T'],tdict['S'],meshes)

            if inargs('jacobian_angle'):
                jacobian_angle.calc(jacobian_depth)

            if inargs('bp'):
                bpd.calc(tdict['ssh'],rho.data)

            if inargs('mont'):
                mgd.calc(tdict['ssh'],tdict['T'],tdict['S'])
                for instance in instancedict.values():
                    instance.calc(mgd)

            for x in args.passive_s:
                instance = passive_s_dict[x]
                trname = x[:-2]
                instance.calc(tdict[trname],mgd)

            if inargs('Sin'):
                Sind.calc(tdict['EmP'],tdict['EmPs'],tdict['sss'])

            if inargs('Win'):
                Wind.calc(tdict['EmP'])

            if inargs('Bin'):
                Bind.calc(tdict['EmPs'],tdict['Hin'],tdict['sst'],tdict['sss'])

            if inargs('THin'):
                THind.calc(tdict['Hin'],tdict['EmP'],tdict['sst'])

            if inargs('IceW'):
                IceWd.calc(tdict['aice'],tdict['hice'])

            if inargs('SnowW'):
                SnowWd.calc(tdict['aice'],tdict['hsnow'])

            for xbar in inargs.arguments:
                if xbar.endswith('bar'):
                    x = xbar[:-3]
                    if nemo_dimensions[x] == 3:
                        bar3d[x].calc(idict[x],tdict['ssh'])
                        print('%s= ' % xbar,bar3d[x].data.nos)
                    elif nemo_dimensions[x] == 2:
                        bar2d[x].calc(idict[x])

            for zm in inargs.arguments:
                if zm.startswith('zm'):
                    x = zm[2:]
                    zm3d[x].calc(idict[x],tdict.get('ssh'))

            for gs in inargs.arguments:
                if gs.startswith('gs'):
                    x = gs[2:]
                    if len(idict[x].nos.shape)==2:
                        gsdict[x].calc(idict[x])
                    elif len(idict[x].nos.shape)==3:
                        gsdict[x].calc(idict[x],tdict['ssh'])

            if inargs('heatsum'):
                G3Heatd.calc(idict['gsT'])

            if inargs('FWsum'):
                G3FWd.calc(idict['gsS'],idict['gsssh'])


    def do_date(p,y,month,day,first):
        if first:
            first = False
        else:
            t01 = time.time()
            for fext in fexts:
                pathname = model_run(passno=p,years=y,month=month,
                                     fext=fext,day5=day)
                cdf_file = DCDF4(pathname)
                for t in tdict.keys():
                    if t in fexttracerd[fext]:
                        cdf_file.update(idict[t])
            do_on_file()
            t03 = time.time()
            print (f'time  to get data from {pathname} is {t03-t01}')
            print (f'time  after  got data from {pathname} is {t03-t00}')

        t01 = time.time()
        print (f'time  before writing {year, month, day} is {t01-t00}')
        threedcdf(xdict,year,month=month,day=day)
        t03 = time.time()
        print (f'time  to perform  write into output file is {t03-t01}')
        print (f'time  after  write into output file is {t03-t00}')


    def do_infile(infile):
        t01 = time.time()
        print (f'time  before processing {infile} etc is {t01-t00}')
        for fext in fexts:

            pathname = infile[:-3]
            while pathname[-1].isalpha():
                pathname = pathname[:-1]

            pathname += fext + '.nc'

            cdf_file = DCDF4(pathname)
            for t in tdict.keys():
                if t in fexttracerd[fext]:
                    cdf_file.update(idict[t])

        do_on_file()
        t02 = time.time()
        print (f'time  to process {pathname} is  {t02-t01}')
        print (f'time  before writing to output file & after'
               f'processing {pathname} is {t02-t00}')
        threedcdf(xdict,None,month=None,day=None)
        t03 = time.time()
        print (f'time  to perform  write into output file is {t03-t02}')
        print (f'time  after  write into output file is {t03-t00}')


    if not use_infile:
        first = True
        for year in years:
            if args.passno is not None:
                dy = year - args.year00
                y, p = dy % 60 + args.year00, dy/60 + passno0
            else:
                y = year
                p = None
            print(year,y, first)
            for month in months:
                print(y,p,month, first)
                do_date(p,y,month,None,first)
                print(y,p,month, first)

            for day in days:
                print(y,p,day, first)
                # explicitly don't call to avoid double call when both month and day are None
                if day is not None:
                    do_date(p,y,None,day,first)
                    print(y,p,day, first)
    else: # infile is not None
        threedcdf(xdict,None,month=None,day=None)
        while args.infile:
            infile = args.infile.pop(0)
            do_infile(infile)
 

    threedcdf.close()
