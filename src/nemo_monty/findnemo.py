# -*- coding: utf-8 -*-
from __future__ import print_function
from os.path import join as pjoin, exists, expanduser, isdir, splitext
import os
import sys
import platform
import time
from argparse import ArgumentParser

def recursive_glob(rootdir='.', suffix='',excludes=[]):
    paths = []
    # domain_ends = ['mask.nc','mesh_zgr.nc','mesh_hgr.nc']
    for rootdir, dirnames, filenames in os.walk(rootdir,followlinks=True):
       for exclude in excludes:
            if exclude in dirnames:
                dirnames.remove(exclude)
       for filename in filenames:
           if filename.endswith(suffix) and not splitext(filename)[0][-4:].isdigit():
               path = pjoin(rootdir, filename)
               paths.append(path)
    return paths

def simple_glob(rootdir='.', suffix='',diagfile=None):
    if diagfile is None:
        print('searching ',rootdir)
    else:
        print ('searching ',rootdir, file = diagfile)

    # print os.listdir(rootdir)
    # paths = [pjoin(rootdir,f) for f in os.listdir(rootdir) if os.path.isfile(pjoin(rootdir,f))]
    paths = [pjoin(rootdir,f) for f in os.listdir(rootdir) if f.endswith(suffix) and
             os.path.isfile(pjoin(rootdir,f)) and not splitext(f)[0][-4:].isdigit()]
    return paths

def get_dirs(basedir):
    if not isdir(basedir):
        return None
    print(basedir)
    print (os.listdir(basedir) )
    return [name for name in os.listdir(basedir)
            if isdir(pjoin(basedir, name))]

class Find_Run(object):
    first = True
    def __init__(self):
        for s in ['centos','macOS']:
            if s in platform.platform():
                system = s
                if Find_Run.first:
                    print('system is ',system)
                break
        else:
            sys.exit('system is %s, not one of %s' % (platform.platform(),' '.join(s)))

        rootdirs = {'centos':['/hpcdata/scratch/omfman/ALTIX1',
                              '/hpcdata/scratch2'],
                    'macOS':['/Volumes/Achernar/Data/NEMO/ORCA12',
                             '/Users/agn/Data/NEMO/ORCA12']}
        self.rootdirs = rootdirs[s]
        if Find_Run.first:
            Find_Run.first = False
            print('rootdirs are %s' % ' '.join(self.rootdirs))

    def find_run(self,runname, rootdirs=None, branchdirs=['NEW_ORCA1']):
        if rootdirs is None:
            rootdirs = self.rootdirs

        for rootdir in rootdirs:
            dirs = get_dirs(rootdir)
            if dirs is None: continue
            if runname in dirs:
                return pjoin(rootdir,runname)
            for branchdir in branchdirs:
                if branchdir in dirs:
                    subdirs = get_dirs(pjoin(rootdir,branchdir))
                    if runname in subdirs:
                        return pjoin(rootdir,branchdir,runname)

        else:
            sys.exit('%s directory not found' % runname)


class Run(object):
    first = True
    def __init__(self, runname, refresh=False, writecache=True, excludes=[], subdir=None):
        self.ignore_pass = False
        excludes += ['WORK','old','WORK_ORCA2','analysis','regions','Hector']
        cachedirs = ['~/.nemocache','~agn/.nemocache']
        for cdir in cachedirs:
            d = expanduser(cdir)
            if exists(d):
                cachedir = d
                break
        else:
            cachedir = expanduser(cachedirs[0])
            os.makedirs(cachedir)

        if subdir is None:
            cachename = runname
        else:
            cachename = '%s+%s' % (runname,subdir)
        cachepath = pjoin(cachedir,cachename)

        cache = refresh or not exists(cachepath)

        find_run = Find_Run()

        if cache:
            rootdir = find_run.find_run(runname)
            if subdir is not None:
                rootdir = pjoin(rootdir,subdir)
            self.valfiles = recursive_glob(rootdir=rootdir, suffix='nc',
                                           excludes=excludes)
            print('found files')
            if writecache:
                with open(cachepath,'w') as f:
                    f.writelines("%s\n" % item for item in self.valfiles )
                print('%s written out' % cachepath)
        else:
            with open(cachepath,'r') as f:
                t0 = time.time()
                self.valfiles = f.read().splitlines()
                t1 = time.time()
                if Run.first:
                    Run.first = False
                    print ('%s read in, took %g\n' % (cachepath, t1-t0))

    def findfiles(self,strings,suffix,nostrings=[]):
        flist =[]
        # print 'In findfiles suffix=', suffix
        # print self.valfiles
        # print [x for x in self.valfiles if x.endswith(suffix)]
        for vf in [x for x in self.valfiles if x.endswith(suffix)]:
            for s in strings:
                if s not in vf:
                    break
            else:
                for s in nostrings:
                    if s in vf:
                        break
                else:
                    flist.append(vf)
        # print flist
        flist = [x for x in flist if not x[x.find(suffix)-1].isalpha()]
        return flist

    def findfile(self,years=None,month=None,day5=None,day1=None,
                 order=['v3.3.1_dfs5.1.1','v3.3.1'], passno=None,fext=''):
        strings = []
        nostrings = []
        if years is not None:
            if not isinstance(years,(list,tuple)):
                years = [years]
            if len(years)>1:
                ystring_dir = '%4ito%4i' % tuple(years)
            else:
                ystring_dir = '%4i' % years[0]
                nostrings.append('to%s' % ystring_dir)

            if 'restart' not in fext:
                ystring = '_%s' % ystring_dir
                if month is None and day5 is None and day1 is None:
                    ystring = '%sy' % ystring_dir

            strings.append(ystring)
        else:
            ystring_dir = ''

        if month is not None:
            strings.append('m%02i' % month)
        if day5 is not None:
            strings.append('%04id05' % day5)
        if day1 is not None:
            strings.append('%04id01' % day1)
        if passno is not None and not self.ignore_pass:
            if 'restart' in fext:
                 if years is not None and len(years)==1:
                     strings.append('P%i_%s' % (passno,ystring))
            else:
                strings.append('PASS%i' % passno)

        suffix = '%s.nc' % fext
        flist = self.findfiles(strings,suffix,nostrings)
        flist2 = flist[:]
        y_order = ['%s_%s' % (ystring_dir, o) for o in order]
        for f in flist:
            dirname,fname = os.path.split(f)
            for n,o in enumerate(order):
                if o in dirname:
                    print (n,o,order,order[n:],order[n:] + [ystring_dir])
                    for o0 in y_order[n+1:] + [ystring_dir]:
                        f0 = pjoin(dirname.replace(y_order[n],o0), fname)
                        if f0 in flist2:
                            flist2.remove(f0)
                    break
        # print flist
        # if len(flist2)<1:
        #     sys.exit('cannot find file from\n%s'% '\n'.join(flist))

        return flist2


class ExtraDir(Run):
    def __init__(self,rootdir,diagfile=None):

        self.valfiles = simple_glob(rootdir=rootdir, suffix='nc',diagfile=diagfile)
        if len(self.valfiles) > 0:
            if diagfile is  None:
                print('found files')
            else:
                print( 'found files', file = diagfile)
        self.ignore_pass = True


class MultiDir():
    def __init__(self,runname,xroots=None,refresh=False,diagfile=None):
        self.runname = runname
        self.xroots = xroots
        self.refresh = refresh
        self.diagfile = diagfile

    def __call__(self,**findargs):
        foundfile = False
        if self.xroots is not None:
            for xroot in self.xroots:
                extra = ExtraDir(xroot,diagfile=self.diagfile)
                plist = extra.findfile(**findargs)
                if len(plist)>0:
                    foundfile = True
                    break

        if not foundfile:
            nemorun = Run(self.runname,refresh=self.refresh)
            plist = nemorun.findfile(**findargs)

        if len(plist)<1:
            nemorun = Run(self.runname,refresh=self.refresh, subdir='analysis')
            plist = nemorun.findfile(**findargs)

        if len(plist)<1:
            print(findargs)
            sys.exit('no file found for the above')

        elif len(plist)>1:
            for p in plist:
                dirname,fname = os.path.split(p)
                for end in ['mask.nc','mesh_zgr.nc','mesh_hgr.nc']:
                    if fname.endswith(end) and len(fname)>len(end):
                        plist.remove(p)
            if len(plist)>1:
                for p in plist:
                    dirname,fname = os.path.split(p)
                    if fname in ['mask.nc','mesh_zgr.nc','mesh_hgr.nc']:
                        break
                else:
                    sys.exit('multiple files found %s' % ' '.join(plist))

        return plist[0]


if __name__ == '__main__':
    parser = ArgumentParser(description='Find NEMO output files')
    #parser.add_argument(dest='file_details',help='Details of file being searched for')
    parser.add_argument('-y',dest='years',help='Year(s)',type=int,nargs= '*',default=None)
    parser.add_argument('-m',dest='month',help='Month',type=int,default=None)
    parser.add_argument('--pass',dest='passno',type=int,help='pass number',default=None)
    parser.add_argument('--day5',dest='day5',help='last 5day of 5-day period',type=int,default=None)
    parser.add_argument('--day1',dest='day1',help='day #',type=int,default=None)
    parser.add_argument('-r',dest='runname',help='Run name',default='ORCA1-N403')
    parser.add_argument('-e',dest='fext',help='File identifier',default='T')
    parser.add_argument('--xroot',dest='xroots',help='Extra root directories',nargs= '*',default=None)
    parser.add_argument('--refresh',dest='refresh',help='refresh cache?',
                        action='store_true',default=False)
    args = parser.parse_args()

    findargs = dict(years = args.years, month = args.month, day5 = args.day5,
                            day1 = args.day1, passno = args.passno,fext = args.fext)

    extradir = MultiDir(args.runname,xroots=args.xroots,refresh=args.refresh)
    filename = extradir(**findargs)
    print(filename)

