# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 23:04:33 2019

@author: lindseykitchell

Series of python functions for digitizing fixed and sliding 3D landmarks 
on .vtk, .ply, or .stl files. 
"""

def organizeSubjsFiles(directory):
    import glob
    import os
    subjlist = []
    filelist=[]
    for file in glob.glob(directory+"/*"):
        subjlist.append(os.path.basename(file)[0:5])
        filelist.append(os.path.abspath(file))
    return subjlist,filelist

def choose_fixedLandmarks(file_name,template=[]):
    import vtk
    from dipy.viz import window
    import numpy as np

    

    # import the 3D surface from file
    if file_name[-3:] == 'vtk':
        object = vtk.vtkPolyDataReader()
    if file_name[-3:] == 'ply':
        object = vtk.vtkPLYReader()
    if file_name[-3:] == 'stl':
        object = vtk.vtkSTLReader()
    object.SetFileName(file_name)
    
    objectMapper = vtk.vtkPolyDataMapper()
    objectMapper.SetInputConnection(object.GetOutputPort())
    objectMapper.ScalarVisibilityOff()
    
    objectActor=vtk.vtkActor()
    objectActor.SetMapper(objectMapper)
    objectActor.GetProperty().SetColor(0.5,0.5,0.5) # grey
    #left over from other work, don't delete yet    
    #objectActor.GetProperty().SetColor(.24, .70, .44) #mediumseagreen
    #objectActor.GetProperty().SetColor(0.498039, 1, 0.831373) #springgreen
    #objectActor.GetProperty().SetColor(color[0],color[1],color[2])
    
    # Attach to a renderer
    #    ren = vtk.vtkRenderer()
    #    ren.AddActor(objectActor)
    #    ren.SetBackground(0.1, 0.1, 0.1)
    
    # Attach to a window
    #    renWin = vtk.vtkRenderWindow()
    #    renWin.AddRenderer(ren)
    #    #renWin.SetWindowName("surface")
    #    renWin.SetSize(500,500)
    
    # Attach to an interactor
    #    iren = vtk.vtkRenderWindowInteractor()
    #    iren.SetRenderWindow(renWin)
    #    style = vtk.vtkInteractorStyleSwitch()
    #    style.SetCurrentStyleToTrackballCamera()
    #    iren.SetInteractorStyle(style)
    #    iren.Initialize()
    #    iren.Start()
    #    
    #   
    
    #create render window    
    renderer = window.Renderer()
    #add 3D surface to window
    renderer.add(objectActor)
    
    #initialize landmarks
    landmarks = []
    
    #initialize dictionary for tying landmarks to annotations
    actortextdict = {}
    
    def mark(x,y,z):
        #function for placing a landmark sphere
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1)
        res = 20
        sphere.SetThetaResolution(res)
        sphere.SetPhiResolution(res)
        sphere.SetCenter(x,y,z)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
    
        marker = vtk.vtkActor()
        marker.SetMapper(mapper)
        renderer.AddActor(marker)
        marker.GetProperty().SetColor( (1,0,0) )
        
        #annotate the mark        
        atext = vtk.vtkVectorText()
        atext.SetText(str(len(landmarks)+1))
        textMapper = vtk.vtkPolyDataMapper()
        textMapper.SetInputConnection(atext.GetOutputPort())
        textActor = vtk.vtkFollower()
        textActor.SetMapper(textMapper)
        textActor.SetScale(3, 3, 3)
        textActor.AddPosition(x,y,z)
        textActor.SetCamera(renderer.GetActiveCamera())
        actortextdict[marker]=textActor
        renderer.AddActor(textActor)
        
        show_m.iren.Render()
    
    def pick_cell(renwinInteractor, event):
        #function for placing a landmark at mouse click position
        x, y = renwinInteractor.GetEventPosition()
    
        picker = vtk.vtkCellPicker()
        picker.PickFromListOn()
        picker.AddPickList(objectActor)
        picker.SetTolerance(0.01)
        picker.Pick(x, y, 0, renderer)
        points = picker.GetPickedPositions()
        numPoints = points.GetNumberOfPoints()
        if numPoints<1: return
        pnt = points.GetPoint(0)
        mark(*pnt)
        landmarks.append(pnt)
    
        
    def remove_cell(renwinInteractor, event):
        #function for removing a landmark at mouse click position
        x, y = renwinInteractor.GetEventPosition()
        
        picker = vtk.vtkPropPicker()
    #    picker.PickFromListOn()
    #    picker.AddPickList(objectActor)
    #    picker.SetTolerance(0.01)
        picker.Pick(x, y, 0, renderer)
        pickedActor = picker.GetActor()
        
        if pickedActor:
            if pickedActor != objectActor:
                mapper = pickedActor.GetMapper()
                inputs = mapper.GetInput()
                point = inputs.GetCenter()
                if point not in template:
                    renderer.RemoveActor(pickedActor)
                    renderer.RemoveActor(actortextdict[pickedActor])
                    show_m.iren.Render()
                    if point in landmarks:
                        landmarks.remove(point)
                
    
    def marktemplate(x,y,z,landnum):
        #function for placing template landmarks
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1)
        res = 20
        sphere.SetThetaResolution(res)
        sphere.SetPhiResolution(res)
        sphere.SetCenter(x,y,z)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
    
        marker = vtk.vtkActor()
        marker.SetMapper(mapper)
        renderer.AddActor(marker)
        marker.GetProperty().SetColor( (1,1,0) )
        
        #annotate the mark
        atext = vtk.vtkVectorText()
        atext.SetText(str(landnum+1))
        textMapper = vtk.vtkPolyDataMapper()
        textMapper.SetInputConnection(atext.GetOutputPort())
        textActor = vtk.vtkFollower()
        textActor.SetMapper(textMapper)
        textActor.SetScale(3, 3, 3)
        textActor.AddPosition(x,y,z)
        textActor.SetCamera(renderer.GetActiveCamera())
        renderer.AddActor(textActor)        
        
        show_m.iren.Render()
    
    
    show_m = window.ShowManager(renderer, size=(800, 800))
    show_m.iren.AddObserver('LeftButtonPressEvent', pick_cell)
    show_m.iren.AddObserver('RightButtonPressEvent', remove_cell)
    
    tempnum = 0
    for lnum in range(len(template)):
        lm = template[lnum]
        marktemplate(lm[0], lm[1], lm[2], tempnum)
        tempnum+=1    
    
    
    show_m.initialize()
    show_m.render()
    show_m.start()
    
    return np.array(landmarks)
    

def show_fixedLandmarks(file_name, landmarks):
    import vtk
    from dipy.viz import window
            
    # Read the surface from file
    if file_name[-3:] == 'vtk':
        object = vtk.vtkPolyDataReader()
    if file_name[-3:] == 'ply':
        object = vtk.vtkPLYReader()
    if file_name[-3:] == 'stl':
        object = vtk.vtkSTLReader()
    object.SetFileName(file_name)
    
    objectMapper = vtk.vtkPolyDataMapper()
    objectMapper.SetInputConnection(object.GetOutputPort())
    objectMapper.ScalarVisibilityOff()
    
    objectActor=vtk.vtkActor()
    objectActor.SetMapper(objectMapper)
    objectActor.GetProperty().SetColor(0.5,0.5,0.5)

    renderer = window.Renderer()
    
    renderer.add(objectActor)

    def mark(x,y,z,landnum):
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1)
        res = 20
        sphere.SetThetaResolution(res)
        sphere.SetPhiResolution(res)
        sphere.SetCenter(x,y,z)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
    
        marker = vtk.vtkActor()
        marker.SetMapper(mapper)
        renderer.AddActor(marker)
        marker.GetProperty().SetColor( (1,0,0) )
        
        #annotate the mark
        atext = vtk.vtkVectorText()
        atext.SetText(str(landnum+1))
        textMapper = vtk.vtkPolyDataMapper()
        textMapper.SetInputConnection(atext.GetOutputPort())
        textActor = vtk.vtkFollower()
        textActor.SetMapper(textMapper)
        textActor.SetScale(3, 3, 3)
        textActor.AddPosition(x,y,z)
        textActor.SetCamera(renderer.GetActiveCamera())
        renderer.AddActor(textActor)        
        
        show_m.iren.Render()

    show_m = window.ShowManager(renderer, size=(800, 800))
    
    for lnum in range(len(landmarks)):
        lm = landmarks[lnum]
        mark(lm[0], lm[1], lm[2], lnum)

    show_m.initialize()
    show_m.render()
    show_m.start()

def choose_curveLandmarks(file_name, fixedLandmarks = []):
    import vtk
    from dipy.viz import window
    import numpy as np
    #import os
    
    
    # Read the surface from file
    if file_name[-3:] == 'vtk':
        object = vtk.vtkPolyDataReader()
    if file_name[-3:] == 'ply':
        object = vtk.vtkPLYReader()
    if file_name[-3:] == 'stl':
        object = vtk.vtkSTLReader()
    object.SetFileName(file_name)
    
    objectMapper = vtk.vtkPolyDataMapper()
    objectMapper.SetInputConnection(object.GetOutputPort())
    objectMapper.ScalarVisibilityOff()
    
    objectActor=vtk.vtkActor()
    objectActor.SetMapper(objectMapper)
    objectActor.GetProperty().SetColor(0.5,0.5,0.5) # grey
    
    renderer = window.Renderer()
    
    renderer.add(objectActor)
    
    landmarks = []
    
    actortextdict = {}
    actorlinedict = {}
    
    
    def mark(x,y,z):
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1)
        res = 20
        sphere.SetThetaResolution(res)
        sphere.SetPhiResolution(res)
        sphere.SetCenter(x,y,z)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
    
        marker = vtk.vtkActor()
        marker.SetMapper(mapper)
        renderer.AddActor(marker)
        marker.GetProperty().SetColor( (0,0,1) )
        
        #annotate the mark        
        atext = vtk.vtkVectorText()
        atext.SetText(str(len(landmarks)+1))
        textMapper = vtk.vtkPolyDataMapper()
        textMapper.SetInputConnection(atext.GetOutputPort())
        textActor = vtk.vtkFollower()
        textActor.SetMapper(textMapper)
        textActor.SetScale(3, 3, 3)
        textActor.AddPosition(x,y,z)
        textActor.SetCamera(renderer.GetActiveCamera())
        actortextdict[marker]=textActor
        renderer.AddActor(textActor)
        
        
        #add line
        if len(landmarks) > 0:
            lineSource = vtk.vtkLineSource()
            lineSource.SetPoint1(landmarks[-1])
            lineSource.SetPoint2([x,y,z])
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(lineSource.GetOutputPort())
            lineActor = vtk.vtkActor()
            lineActor.SetMapper(mapper)
            lineActor.GetProperty().SetLineWidth(4)
            lineActor.GetProperty().SetColor((0,1,0))
            actorlinedict[marker]=lineActor
            renderer.AddActor(lineActor)
        
        show_m.iren.Render()
    
    def pick_cell(renwinInteractor, event):
    
        x, y = renwinInteractor.GetEventPosition()
    
        picker = vtk.vtkCellPicker()
        picker.PickFromListOn()
        picker.AddPickList(objectActor)
        picker.SetTolerance(0.01)
        picker.Pick(x, y, 0, renderer)
        points = picker.GetPickedPositions()
        numPoints = points.GetNumberOfPoints()
        if numPoints<1: return
        pnt = points.GetPoint(0)
        mark(*pnt)
        landmarks.append(pnt)
    
        
    def remove_cell(renwinInteractor, event):
        x, y = renwinInteractor.GetEventPosition()
        
        picker = vtk.vtkPropPicker()
    #    picker.PickFromListOn()
    #    picker.AddPickList(objectActor)
    #    picker.SetTolerance(0.01)
        picker.Pick(x, y, 0, renderer)
        pickedActor = picker.GetActor()
        
        if pickedActor:
            if pickedActor != objectActor:
                mapper = pickedActor.GetMapper()
                inputs = mapper.GetInput()
                point = inputs.GetCenter()
                if point not in fixedLandmarks:
                    renderer.RemoveActor(pickedActor)
                    renderer.RemoveActor(actortextdict[pickedActor])
                    if len(landmarks) > 1:
                        renderer.RemoveActor(actorlinedict[pickedActor])
                    show_m.iren.Render()
                    if point in landmarks:
                        landmarks.remove(point)
            
    def markfixed(x,y,z,landnum):
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1)
        res = 20
        sphere.SetThetaResolution(res)
        sphere.SetPhiResolution(res)
        sphere.SetCenter(x,y,z)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
    
        marker = vtk.vtkActor()
        marker.SetMapper(mapper)
        renderer.AddActor(marker)
        marker.GetProperty().SetColor( (1,0,0) )
        
        #annotate the mark
        atext = vtk.vtkVectorText()
        atext.SetText(str(landnum+1))
        textMapper = vtk.vtkPolyDataMapper()
        textMapper.SetInputConnection(atext.GetOutputPort())
        textActor = vtk.vtkFollower()
        textActor.SetMapper(textMapper)
        textActor.SetScale(3, 3, 3)
        textActor.AddPosition(x,y,z)
        textActor.SetCamera(renderer.GetActiveCamera())
        renderer.AddActor(textActor)        
        
        show_m.iren.Render()
    

    
    show_m = window.ShowManager(renderer, size=(800, 800))
    show_m.iren.AddObserver('LeftButtonPressEvent', pick_cell)
    show_m.iren.AddObserver('RightButtonPressEvent', remove_cell)
    
    #    if fixedLandmarks:
    for lnum in range(len(fixedLandmarks)):
        lm = fixedLandmarks[lnum]
        markfixed(lm[0], lm[1], lm[2], lnum)
        
    show_m.initialize()
    show_m.render()
    show_m.start()
    
    return np.array(landmarks)
    

def show_curveLandmarks(file_name, landmarks):
    import vtk
#    from vtk.util.numpy_support import vtk_to_numpy
    from dipy.viz import window
    #import os
    
    
    
    # Read the surface from file
    if file_name[-3:] == 'vtk':
        object = vtk.vtkPolyDataReader()
    if file_name[-3:] == 'ply':
        object = vtk.vtkPLYReader()
    if file_name[-3:] == 'stl':
        object = vtk.vtkSTLReader()
    object.SetFileName(file_name)
    
    objectMapper = vtk.vtkPolyDataMapper()
    objectMapper.SetInputConnection(object.GetOutputPort())
    objectMapper.ScalarVisibilityOff()
    
    objectActor=vtk.vtkActor()
    objectActor.SetMapper(objectMapper)
    objectActor.GetProperty().SetColor(0.5,0.5,0.5)

    renderer = window.Renderer()
    
    renderer.add(objectActor)

    def mark(x,y,z,landnum):
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1)
        res = 20
        sphere.SetThetaResolution(res)
        sphere.SetPhiResolution(res)
        sphere.SetCenter(x,y,z)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
    
        marker = vtk.vtkActor()
        marker.SetMapper(mapper)
        renderer.AddActor(marker)
        marker.GetProperty().SetColor( (0,0,1) )
        
        #annotate the mark
        atext = vtk.vtkVectorText()
        atext.SetText(str(landnum+1))
        textMapper = vtk.vtkPolyDataMapper()
        textMapper.SetInputConnection(atext.GetOutputPort())
        textActor = vtk.vtkFollower()
        textActor.SetMapper(textMapper)
        textActor.SetScale(3, 3, 3)
        textActor.AddPosition(x,y,z)
        textActor.SetCamera(renderer.GetActiveCamera())
        renderer.AddActor(textActor)        
        
        #add line
        if landnum > 0:
            lineSource = vtk.vtkLineSource()
            lineSource.SetPoint1(landmarks[landnum-1])
            lineSource.SetPoint2([x,y,z])
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(lineSource.GetOutputPort())
            lineActor = vtk.vtkActor()
            lineActor.SetMapper(mapper)
            lineActor.GetProperty().SetLineWidth(4)
            lineActor.GetProperty().SetColor((0,1,0))
            renderer.AddActor(lineActor)
            
        show_m.iren.Render()

    show_m = window.ShowManager(renderer, size=(800, 800))
    
    for lnum in range(len(landmarks)):
        lm = landmarks[lnum]
        mark(lm[0], lm[1], lm[2], lnum)

    show_m.initialize()
    show_m.render()
    show_m.start()


def show_allLandmarks(file_name, fixedLandmarks, curveLandmarks):
    # input the curve landmarks in a list: curveLandmarks = [curve1, curve2]    
    #if you dont wish to input both fixed and curve landmarks you can 
    # place a [] as the input
    import vtk
#    from vtk.util.numpy_support import vtk_to_numpy
    from dipy.viz import window
    #import os
    
    
    
    # Read the surface from file
    if file_name[-3:] == 'vtk':
        object = vtk.vtkPolyDataReader()
    if file_name[-3:] == 'ply':
        object = vtk.vtkPLYReader()
    if file_name[-3:] == 'stl':
        object = vtk.vtkSTLReader()
    object.SetFileName(file_name)
    
    objectMapper = vtk.vtkPolyDataMapper()
    objectMapper.SetInputConnection(object.GetOutputPort())
    objectMapper.ScalarVisibilityOff()
    
    objectActor=vtk.vtkActor()
    objectActor.SetMapper(objectMapper)
    objectActor.GetProperty().SetColor(0.5,0.5,0.5)

    renderer = window.Renderer()
    
    renderer.add(objectActor)

    def markcurve(x,y,z,landnum, curve, annotnum):
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1)
        res = 20
        sphere.SetThetaResolution(res)
        sphere.SetPhiResolution(res)
        sphere.SetCenter(x,y,z)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
    
        marker = vtk.vtkActor()
        marker.SetMapper(mapper)
        renderer.AddActor(marker)
        marker.GetProperty().SetColor( (0,0,1) )
        
        #annotate the mark
        atext = vtk.vtkVectorText()
        atext.SetText(str(annotnum+1))
        textMapper = vtk.vtkPolyDataMapper()
        textMapper.SetInputConnection(atext.GetOutputPort())
        textActor = vtk.vtkFollower()
        textActor.SetMapper(textMapper)
        textActor.SetScale(3, 3, 3)
        textActor.AddPosition(x,y,z)
        textActor.SetCamera(renderer.GetActiveCamera())
        renderer.AddActor(textActor)        
        
        #add line
        if landnum > 0:
            lineSource = vtk.vtkLineSource()
            lineSource.SetPoint1(curve[landnum-1])
            lineSource.SetPoint2([x,y,z])
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(lineSource.GetOutputPort())
            lineActor = vtk.vtkActor()
            lineActor.SetMapper(mapper)
            lineActor.GetProperty().SetLineWidth(4)
            lineActor.GetProperty().SetColor((0,1,0))
            renderer.AddActor(lineActor)
            
        show_m.iren.Render()

    def markfixed(x,y,z,landnum):
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1)
        res = 20
        sphere.SetThetaResolution(res)
        sphere.SetPhiResolution(res)
        sphere.SetCenter(x,y,z)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
    
        marker = vtk.vtkActor()
        marker.SetMapper(mapper)
        renderer.AddActor(marker)
        marker.GetProperty().SetColor( (1,0,0) )
        
        #annotate the mark
        atext = vtk.vtkVectorText()
        atext.SetText(str(landnum+1))
        textMapper = vtk.vtkPolyDataMapper()
        textMapper.SetInputConnection(atext.GetOutputPort())
        textActor = vtk.vtkFollower()
        textActor.SetMapper(textMapper)
        textActor.SetScale(3, 3, 3)
        textActor.AddPosition(x,y,z)
        textActor.SetCamera(renderer.GetActiveCamera())
        renderer.AddActor(textActor)        
        
        show_m.iren.Render()

    show_m = window.ShowManager(renderer, size=(800, 800))
    
    totalnum = 0
#    if fixedLandmarks:
    for lnum in range(len(fixedLandmarks)):
        lm = fixedLandmarks[lnum]
        markfixed(lm[0], lm[1], lm[2], totalnum)
        totalnum+=1
    for c in curveLandmarks:
        for lnum in range(len(c)):
            lm = c[lnum]
            markcurve(lm[0], lm[1], lm[2], lnum, c, totalnum)
            totalnum +=1
            
    show_m.initialize()
    show_m.render()
    show_m.start()
    
def resampleCurve(curve, npoints, includeFixed=False, fixedLandmarks = [], fixedStart=1, fixedEnd=1):
    #returns the curve resampled to npoints number of equidistant points
    #if fixed landmarks make up the end points of a curve you can indicate 
    #them with the optional arguments
    #the returned resampled curve will not include the fixed landmark endpoints
    import numpy as np
    from dipy.tracking.streamline import set_number_of_points
    if includeFixed:
        curve = np.vstack((fixedLandmarks[fixedStart-1], curve, fixedLandmarks[fixedEnd-1]))
    resampledcurve = set_number_of_points(curve, npoints)
    if includeFixed:
        resampledcurve=resampledcurve[1:-1]
    return resampledcurve
    
def saveLandmarks_tps(tpsfilename,subjlist,fixedLandmarksdict, curveLandmarksdict):
    tpsfile = open(tpsfilename+'.tps','w') 
    for subj in range(len(subjlist)):
        #get total number of landmarks        
        nlandmarks = len(fixedLandmarksdict[subjlist[subj]])        
        for c in curveLandmarksdict[subjlist[subj]]:
            nlandmarks = nlandmarks + len(c)
        
        tpsfile.write('LM3='+str(nlandmarks)+'\n')
        for i in range(len(fixedLandmarksdict[subjlist[subj]])):
            tpsfile.write(str(fixedLandmarksdict[subjlist[subj]][i][0])+' '+str(fixedLandmarksdict[subjlist[subj]][i][1])+' '+str(fixedLandmarksdict[subjlist[subj]][i][2])+'\n')
        for c in curveLandmarksdict[subjlist[subj]]:
            for i in range(len(c)):
               tpsfile.write(str(c[i][0])+' '+str(c[i][1])+' '+str(c[i][2])+'\n') 
        tpsfile.write('ID='+subjlist[subj]+'\n'+'\n')
    tpsfile.close()

def saveLandmarks_pkl(pklfilename, subjlist, fixedLandmarksDict,curvelandmarksDict):
    import pickle
    pkldict = {}
    pkldict['subjlist']=subjlist
    pkldict['fixedLandmarksDict']=fixedLandmarksDict
    pkldict['curvelandmarksDict']=curvelandmarksDict
    
    pickle.dump( pkldict, open( pklfilename+".p", "wb" ) )

def saveLandmarks_vtk(subjlist, landmarkDict):
    #saves .vtk for single subject
    import vtk 
    
    for i in range(len(subjlist)):
        landmarks = landmarkDict[subjlist[i]]        
        vtkfilename = subjlist[i]+'_landmarks.vtk'        
        
        Points = vtk.vtkPoints()
        for i in range(len(landmarks)):
            Points.InsertNextPoint(landmarks[i][0],landmarks[i][1],landmarks[i][2])
        
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(Points)
        writer = vtk.vtkDataSetWriter()
        writer.SetFileName(vtkfilename)
        writer.SetInputData(polydata)
        writer.Write()

def loadLandmarks_pkl(pklfilename):
    import pickle
    pkldict = pickle.load( open(pklfilename , "rb" ) )
    subjlist = pkldict['subjlist']
    fixedLandmarksDict = pkldict['fixedLandmarksDict']
    curvelandmarksDict = pkldict['curvelandmarksDict']
    return subjlist, fixedLandmarksDict, curvelandmarksDict
    
def loadLandmarks_tps(tpsfilename):
    import numpy as np

    """
    Function to read a .TPS file
    Args:
        tpsfilename (str): path to the .TPS file
    Returns:
        lm (str list): info extracted from 'LM=' field
        im (str list): info extracted from 'IMAGE=' field
        id (str list): info extracted from 'ID=' filed
        coords: returns a 3D numpy array if all the individuals have same
                number of landmarks, otherwise returns a list containing 2d
                matrices of landmarks
    edited from this source: 
    https://gist.github.com/jinyung/1b8fe5735fbfdf07378197cc4c9acc3a
    """

    # open the file
    tps_file = open(tpsfilename, 'r')  # 'r' = read
    tps = tps_file.read().splitlines()  # read as lines and split by new lines
    tps_file.close()

    # initiate lists to take fields of "LM=","IMAGE=", "ID=" and the coords
    lm, im, ID, coords_array = [], [], [], []

    # looping thru the lines
    for i, ln in enumerate(tps):

        # Each individual starts with "LM="
        if ln.startswith("LM"):
            # number of landmarks of this ind
            lm_num = int(ln.split('=')[1])
            # fill the info to the list for all inds
            lm.append(lm_num)
            # initiate a list to take 2d coordinates
            coords_mat = []

            # fill the coords list by reading next lm_num of lines
            for j in range(i + 1, i + 1 + lm_num):
                coords_mat.append(tps[j].split(' '))  # split lines into values

            # change the list into a numpy matrix storing float vals
            coords_mat = np.array(coords_mat, dtype=float)
            # fill the ind 2d matrix into the 3D coords array of all inds
            coords_array.append(coords_mat)
            # coords_array.append(coords_mat)

        # Get info of IMAGE= and ID= fields
        if ln.startswith("IMAGE"):
            im.append(ln.split('=')[1])

        if ln.startswith("ID"):
            ID.append(ln.split('=')[1])

    # check if all inds contains same number of landmarks
#    all_lm_same = all(x == lm[0] for x in lm)
    # if all same change the list into a 3d numpy array
#    if all_lm_same:
#        coords_array = np.dstack(coords_array)

    subjlist = ID
    allLandmarksDict = {}
    for i in range(len(subjlist)):
        allLandmarksDict[subjlist[i]] = coords_array[i]
    # return results in dictionary form
    return subjlist, allLandmarksDict
  
def define_slidingLandmarks(fixedLandmarks, curveLandmarks,fixedStart_list, fixedEnd_list):
    import pandas as pd 
    import numpy as np    
    curves = []

    nlandmarks = len(fixedLandmarks)
    for c in range(len(curveLandmarks)):
        if fixedStart_list[c]:
            nlandmarks+=1
            curves.append([fixedStart_list[c],nlandmarks,nlandmarks+1])
            for i in range(len(curveLandmarks[c])-2):
                nlandmarks+=1
                curves.append([nlandmarks-1,nlandmarks,nlandmarks+1])
            nlandmarks+=1
            curves.append([nlandmarks-1,nlandmarks,fixedEnd_list[c]])   
        else:
            for i in range(len(curveLandmarks[c])-2):
                nlandmarks+=1
                curves.append([nlandmarks,nlandmarks+1,nlandmarks+2])
            nlandmarks+=2
    df = pd.DataFrame(np.array(curves))
    df.to_csv("curveslide.csv",header=['before', 'slide', 'after'])
    return curves
            
  
def plotAllSubjectLandmarks(subjlist,fixedLandmarkDict, curveLandmarkDict, treatfixedascurve=False, sphereradius=1, group=[]):
    import vtk
    from dipy.viz import window
    import numpy as np

    def markcurve(x,y,z,landnum, curve, color):
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(sphereradius)
        res = 20
        sphere.SetThetaResolution(res)
        sphere.SetPhiResolution(res)
        sphere.SetCenter(x,y,z)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
    
        marker = vtk.vtkActor()
        marker.SetMapper(mapper)
        renderer.AddActor(marker)
        marker.GetProperty().SetColor( color )
                
        
        #add line
        if landnum > 0:
            lineSource = vtk.vtkLineSource()
            lineSource.SetPoint1(curve[landnum-1])
            lineSource.SetPoint2([x,y,z])
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(lineSource.GetOutputPort())
            lineActor = vtk.vtkActor()
            lineActor.SetMapper(mapper)
            lineActor.GetProperty().SetLineWidth(4)
            lineActor.GetProperty().SetColor(color)
            renderer.AddActor(lineActor)
            
        show_m.iren.Render()

    def markfixed(x,y,z, color):
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(sphereradius)
        res = 20
        sphere.SetThetaResolution(res)
        sphere.SetPhiResolution(res)
        sphere.SetCenter(x,y,z)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())
    
        marker = vtk.vtkActor()
        marker.SetMapper(mapper)
        renderer.AddActor(marker)
        marker.GetProperty().SetColor( color )
        
        show_m.iren.Render()

    renderer = window.Renderer()
    show_m = window.ShowManager(renderer, size=(800, 800))
    
    colorlist = []
    if group:
        ngroups = len(np.unique(group))
        grpcolors = []
        for i in range(ngroups):
            grpcolors.append(list(np.random.choice(range(256), size=3)/256.))
        for i in range(len(group)):
            colorlist.append(grpcolors[list(np.unique(group)).index(group[i])])
    else:        
        for i in range(len(subjlist)):   
            colorlist.append(list(np.random.choice(range(256), size=3)/256.))    
    
    for subj in range(len(subjlist)):
        
        color = colorlist[subj]

        if type(fixedLandmarkDict) is dict:
            if treatfixedascurve:        
                for lnum in range(len(fixedLandmarkDict[subjlist[subj]])):
                    lm = fixedLandmarkDict[subjlist[subj]][lnum]
                    markcurve(lm[0], lm[1], lm[2], lnum, fixedLandmarkDict[subjlist[subj]], color)
            else:
                for lnum in range(len(fixedLandmarkDict[subjlist[subj]])):
                    lm = fixedLandmarkDict[subjlist[subj]][lnum]
                    markfixed(lm[0], lm[1], lm[2], color)
        if type(curveLandmarkDict) is dict:
            for c in curveLandmarkDict[subjlist[subj]]:
                for lnum in range(len(c)):
                    lm = c[lnum]
                    markcurve(lm[0], lm[1], lm[2], lnum, c, color)

                
    show_m.initialize()
    show_m.render()
    show_m.start()
        
