# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 13:48:38 2019

@author: lindseykitchell
"""

from digitize3DLandmarks import *


subj=9

fixedLandmarkDict[subjlist[subj]]=choose_fixedLandmarks(filelist[subj])

curve1 = choose_curveLandmarks(filelist[subj],fixedLandmarkDict[subjlist[subj]])

curve1res = resampleCurve(curve1,22,includeFixed=True, fixedLandmarks=fixedLandmarkDict[subjlist[subj]],fixedStart=1,fixedEnd=12)

curve2 = choose_curveLandmarks(filelist[subj],fixedLandmarkDict[subjlist[subj]])

curve2res = resampleCurve(curve2,22,includeFixed=True, fixedLandmarks=fixedLandmarkDict[subjlist[subj]],fixedStart=1,fixedEnd=12)

curveLandmarkDict[subjlist[subj]]=[curve1res, curve2res]

show_allLandmarks(filelist[subj],fixedLandmarkDict[subjlist[subj]], curveLandmarkDict[subjlist[subj]])


saveLandmarks_tps('left_cingulumcingulate',subjlist,fixedLandmarkDict, curveLandmarkDict)


define_slidingLandmarks(fixedLandmarkDict[subjlist[0]], curveLandmarkDict[subjlist[0]],[1,1], [12,12])