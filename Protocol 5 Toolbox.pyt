# -*- coding: utf-8 -*-
'''
Developed by Brennan Smith and Joshua Thompson (Anne Arundel County BWPR)
User provides a line feature (use template!) with required p5 fields populated. Creates Equilibrium channel DEM to compare against County DEM.

pwsmit32@aacounty.org
pwthom19@aacounty.org

2022-06-10      v2.1        GitHib initial
2024-06-12      v2.2        Don't allow modeled equilibrium DEM to be higher than existing reference DEM (BS)
'''

import arcpy
from arcpy import env
from arcpy.sa import *
from arcpy.ia import * 
import os
import pandas as pd
import numpy as np
import math
import statistics as stat

class Toolbox(object):
    def __init__(self):
        self.label = "Protocol 5 Toolbox"
        self.alias = "p5t"
        self.tools = [P5channelWeqSlope,setRasStretch]


class P5channelWeqSlope(object):
    def __init__(self):
        self.label = "Protocol 5 Equilibrium Channel w/ Eq Slope Calc"
        self.description = "Returns the Protocol 5 Equilibrium Channel DEM for input line template"
        self.canRunInBackground = True

    def getParameterInfo(self):
        #Define parameter definitions
        params = []
        
        #params[0] 
        in_line = arcpy.Parameter(
            displayName="Input line template",
            name="in_line",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        in_line.filter.list = ["Polyline"]  #must be line feature
        
        #params[1]
        in_raster_ref = arcpy.Parameter(
            displayName="Input reference DEM (ft)",
            name="in_raster_ref",
            datatype=["GPRasterLayer","DERasterDataset"],
            parameterType="Required",
            direction="Input")
        in_raster_ref.value = r"https://gis.aacounty.org/image/services/DEM/DEM_2020/ImageServer"
        
        #params[2]
        start_Z = arcpy.Parameter(
            displayName="Reach connectivity method",
            name="start_Z",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        start_Z.filter.list = ['Segmented (use DEM for initial Z)', 'Continuous (link using fromOID)']
        start_Z.value = 'Segmented (use DEM for initial Z)'
        
        #params[3]
        out_ras = arcpy.Parameter(
            displayName="Output equilibrium raster",
            name="out_ras",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")
        out_ras.value = arcpy.CreateUniqueName("EQraster", arcpy.env.workspace)
        
        #params[4] 
        out_poly = arcpy.Parameter(
            displayName="Output equilibrium polygon",
            name="out_poly",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Output")
        out_poly.value = arcpy.CreateUniqueName("EQpoly", arcpy.env.workspace)
        
        #params[5] 
        flip_check = arcpy.Parameter(
            displayName="Check for / enforce uphill line direction",
            name="flip_check",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        flip_check.value = False      
        
        #params[6] 
        calc_eq = arcpy.Parameter(
            displayName="Recalculate equilibrium slope",
            name="calc_eq",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        calc_eq.value = True   
        
        params =   [in_line,
                    in_raster_ref,
                    start_Z,
                    out_ras,
                    out_poly,
                    flip_check,
                    calc_eq]
        
        return params

    def isLicensed(self):
    # Tool needs 3D and Spatial Analyst
        try:
            if arcpy.CheckExtension("3D") != "Available":
                raise Exception
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False  
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        p = {p.name: p for p in parameters}
           
        if p['in_line'].altered:
            
            in_line =  p['in_line'].value
            reqFields = ['eqBankSlope', 'exBottomWidth', 'exBankHeight','bedType']
            
            # Throw a warning if exBedSlope values are missing
            exBeds = [row[0] for row in arcpy.da.SearchCursor(in_line,'exBedSlope')]
            if not all(exBeds):
                p['in_line'].setWarningMessage('Missing existing bed slope values will be computed from DEM')
                msgFlag=1
                
            # loop through each input row, check bed type, update required field list, and populate message string if missing data
            msg = ''
            with arcpy.da.SearchCursor(in_line, [arcpy.Describe(in_line).OIDFieldName,'bedType']) as cursor:
                for row in cursor:
                    checkFields = reqFields.copy()
                    bedType = row[1]
                    if bedType == 'Custom':
                        checkFields.append('eqBedSlope')
                    elif bedType == 'Cohesive':
                        checkFields.append('drainageAcres')                        
                    elif bedType == 'Sand and Gravel':
                        checkFields.extend(['manningsN', 'flowRate', 'initialValDepth', 'exBankSlope'])
                    elif bedType == 'Coarser Than Sand': 
                        checkFields.extend(['manningsN','shieldsParam', 'critBedParam','chanFormDischperUnitWidth','designDisch','meanGrainSize','medGrainSize'])
                    
                    # now check that all the required fields have values
                    getRow = [r for r in arcpy.da.SearchCursor(in_line,checkFields,where_clause = arcpy.Describe(in_line).OIDFieldName+"="+str(row[0]))][0]
                    if not all(getRow):
                        missingFields = [checkFields[f] for f in [i for i,x in enumerate(getRow) if not x]]
                        msg += 'Input feature {} has zero or missing values for required fields: {} \n'.format(row[0],missingFields)
            
            if len(msg) > 0:
                p['in_line'].setErrorMessage(msg)
                
            # Throw a warning if fromOID is empty and the Continuous option is turned on    
            if p['start_Z'].value == 'Continuous (link using fromOID)':
                oids = [row[0] for row in arcpy.da.SearchCursor(in_line,'fromOID')]
                if not any(oids):
                    p['start_Z'].setWarningMessage('No values found for field [fromOID]. Reach elevation will not be continuous')
        
            # If calc_eq is unchecked, make sure there are eqBedSlope values available.
            if p['calc_eq'].value is False:
                eqBeds = [row[0] for row in arcpy.da.SearchCursor(in_line,'eqBedSlope')]
                if not all(eqBeds):
                    p['calc_eq'].setWarningMessage('Missing equilibrium bed slope values will be calculated anyway')
        
        return

    def execute(self, parameters, messages):
    
        arcpy.env.addOutputsToMap = 0
        arcpy.env.overwriteOutput = True
        
        # get parameters
        p = {p.name: p for p in parameters}
        in_line         = p['in_line'].valueAsText
        in_raster_ref   = p['in_raster_ref'].valueAsText
        start_Z         = p['start_Z'].valueAsText
        out_ras         = p['out_ras'].valueAsText
        out_poly        = p['out_poly'].valueAsText
        flip_check      = p['flip_check'].value
        calc_eq         = p['calc_eq'].value
        
        # To make sure we keep only the selected input features, make a new feature layer
        # Not sure why it wasn't honoring input selections in the first place...
        oids = [row[0] for row in arcpy.da.SearchCursor(in_line,arcpy.Describe(in_line).OIDFieldName)]
        arcpy.management.MakeFeatureLayer(in_line, 'line_lyr', 
                                    where_clause=arcpy.Describe(in_line).OIDFieldName+" IN ({:s})".format(','.join(f"{x}" for x in oids)))
                                    
        ## Define Functions ##
        # get generalized channel geometry
        def channel_geom(flowDepth, chanBotWidth, chanSideSlope):
            # channel geometry and normal depth functions taken from rivr R package
            # original code used H:V side slopes. Take inverse of our V:H
            chanSideSlope = 1/chanSideSlope
            A = flowDepth*chanBotWidth + chanSideSlope*flowDepth*flowDepth #A 
            P = chanBotWidth + 2*flowDepth*math.sqrt(chanSideSlope*chanSideSlope + 1) #P
            R = A/P #R
            dAdy = chanBotWidth + 2*flowDepth*chanSideSlope #dAdy
            dPdy = 2*math.sqrt(chanSideSlope*chanSideSlope + 1) #dPdy  
            dTdy = 2*chanSideSlope #dTdy
            dRdy = dAdy/P - A*dPdy/(P*P) #dRdy
            DH = A/dAdy #DH
            ybar = flowDepth*(2*chanBotWidth + dAdy)/(3*(chanBotWidth + dAdy)) # ybar
            chan_geom = {'A': [A], 'P': [P], 'R':[R], 'dAdy': [dAdy], 'dTdy': [dTdy], 'dPdy': [dPdy], 'dRdy': [dRdy], 'DH': [DH], 'ybar': [ybar]}
            return chan_geom
        
        # get normal depth 
        def normal_depth(chanSlope, manningsN, flowRate, initialValDepth, chanBotWidth, chanSideSlope):
            tol = 0.00001
            maxit = 1000 # solve for like 1000 times
            dy = 9999
            i = 0
            while abs(dy) > tol and i < maxit:
                cgm = channel_geom(initialValDepth, chanBotWidth, chanSideSlope)
                dy = (pow(sum(cgm['A']), 5.0/3.0)/pow(sum(cgm['P']), 2.0/3.0) - manningsN*flowRate/(1.486*math.sqrt(chanSlope)))/(sum(cgm['dAdy'])*(5.0/3.0)*pow(sum(cgm['A'])/sum(cgm['P']), 2.0/3.0)-(2.0/3.0)*sum(cgm['dPdy'])*pow(sum(cgm['A'])/sum(cgm['P']), 5.0/3.0))
                initialValDepth = initialValDepth - dy
                i+1
            return initialValDepth
        
        # if bedtype == "Cohesive" calculate eq. bed slope
        def cohesive_EqBedSlope(DA_ac): 
            return 0.0028*(DA_ac/247.105)**-0.33
        
        # if bedtype == "Sand and Gravel" calculate eq. bed slope
        def sand_gravel_EqBedSlope(normDepth): 
            return 0.06/(normDepth*62.43)
        
        # if bedtype == "Coarser Than Sand" calculate eq. bed slope
        def coarserThanSand_EqBedSlope(shieldsParam, critBedParam, chanFormDischperUnitWidth, manningsN, meanGrainSize, designDisch, medGrainSize): 
            return stat.mean([((shieldsParam*critBedParam*1.65)**(10/7))*((1.486/(chanFormDischperUnitWidth*manningsN))**(6/7)), 60.1*((medGrainSize**(10/7))*(manningsN**(9/7)))/((critBedParam**(5/14))*(chanFormDischperUnitWidth**(6/7))), 0.44*(designDisch**-0.46)*(medGrainSize**1.15)])

        # express angle between two points as compass bearing
        def getBearing(x1,y1,x2,y2):
            dX = x2-x1
            dY = y2-y1
            angle = math.atan2(dX, dY)/math.pi*180
            if angle < 0:
                angle = angle+360
            if angle > 360:
                angle = angle - 360
            return angle
        
        # get the difference between two angles/bearings in degrees, expressed within +- 180
        def angleDiff(a1,a2):
            a1 = math.radians(a1)
            a2 = math.radians(a2)
            aDiff = math.atan2(math.sin(a1-a2), math.cos(a1-a2))
            return math.degrees(aDiff)
            
        # rotate a set of XY points: from stack overflow https://stackoverflow.com/a/58781388/6140341
        def rotate(xy, origin=(0, 0), degrees=0):
            ang = np.deg2rad(degrees)
            R = np.array([[np.cos(ang), -np.sin(ang)],
                          [np.sin(ang),  np.cos(ang)]])
            o = np.atleast_2d(origin)
            y = np.atleast_2d(xy)
            return np.squeeze((R @ (xy.T-o.T) + o.T).T)
            
        # Construct a point cloud (array) of the trapezoidal equilibrium channel
        # Using channel parameters, a starting elevation, the XY start/finish points, and some angles: 
        def ReturnPointCloudArray(EqBedSlope,EqBankSlope,ExBottomWidth,ExBankHeight,stX,stY,enX,enY,stZ,exZadj,slopeDiff,angle,lastAngle):
            
            #slopes to angles
            EqBankAngle = math.atan(EqBankSlope) * 180 / math.pi
            
            #get segment length
            dX = enX - stX
            dY = enY - stY
            Length = math.sqrt(dX**2 + dY**2)
            
            ## Build the trapezoid: x is channel width, y is length distance, z is height
            #all vectors are in 1 feet increments, plus ending points
            Xlist = []
            Ylist = []
            Zlist = []
            
            #get the vector of Y and Z values, moving upstream from channel bed center (x=0)
            #these are our longitudinal profile center points
            lp_yvals = list(range(0,math.ceil(Length)))
            lp_yvals.extend([Length])
            lp_zvals = list(np.interp(lp_yvals,[0,Length],[0,Length*EqBedSlope]))
            lp_xvals = [0] * len(lp_yvals)
            #add to master list
            Xlist.extend(lp_xvals)
            Ylist.extend(lp_yvals)
            Zlist.extend(lp_zvals)
            
            ##loop through each y,z long pro pair and build a cross section.
            # Channel Bottom doesn't change so we can define it outside the loop.
            # 1 to the last whole number (even though it's ceil, python range stops one short)
            bottom_xvals = list(range(1,math.ceil(ExBottomWidth/2)))
            # add in the reverse (left side vs right)
            bottom_xvals.extend([x*-1 for x in bottom_xvals])
            # add in the corners, and a zero center point
            bottom_xvals.extend([ExBottomWidth/2,-ExBottomWidth/2,0])
            bottom_xvals.sort()
            
            for y, z in zip(lp_yvals, lp_zvals):
                #set bottom y and z vectors to long pro y,z
                bottom_yvals = [y] * len(bottom_xvals)
                bottom_zvals = [z] * len(bottom_xvals)
                
                ## Channel Banks 
                # height is a function of exHeight, and difference between Eq and Ex beds
                Height = ExBankHeight + exZadj + (y*slopeDiff)
                TopWidth = ((Height*math.tan((90-EqBankAngle)*math.pi/180))*2)+ExBottomWidth
                # build bank X values
                bank_xvals = list(range(math.ceil(ExBottomWidth/2),math.ceil(TopWidth/2)))
                bank_xvals.extend([x*-1 for x in bank_xvals])
                bank_xvals.extend([TopWidth/2,-TopWidth/2])
                bank_xvals.sort()
                # interpolate Z vals (height) of the bank for the X vals of the bank, using our four corner points
                bank_zvals = list(np.interp(bank_xvals,
                                    [-TopWidth/2,-ExBottomWidth/2,ExBottomWidth/2,TopWidth/2],
                                    [Height+(y*EqBedSlope),(y*EqBedSlope),(y*EqBedSlope),Height+(y*EqBedSlope)]))
                bank_yvals = [y] * len(bank_xvals)
                
                # Add to master XYZ Lists
                Xlist.extend(bottom_xvals)
                Ylist.extend(bottom_yvals)
                Zlist.extend(bottom_zvals)
                Xlist.extend(bank_xvals)
                Ylist.extend(bank_yvals)
                Zlist.extend(bank_zvals)
            
            ## compare angle to lastAngle, and build additional cross sections 
            # in X degree increments at the beginning of the segment so we don't have big gaps in our point clouds at each bend
            aDiff   = angleDiff(angle,lastAngle)
            aDiffInc= 5
            aSign   = np.sign(aDiff)
            # if we need to build these extra XS based on our threshold angle (5 degs)
            if abs(aDiff) >= aDiffInc:
                # define the y=0 XS (same as code block above)
                y = 0
                z = 0
                # except only include X values for half the XS (the outside bend, we don't need additional points for the inside bend)
                # right turns are positive aDiff values, where we'd only want the left bank (negative Xvals), so we use -aSign adjustments throughout
                bottom_xvals = list(range(1,math.ceil(ExBottomWidth/2)))
                bottom_xvals.extend([ExBottomWidth/2,0])
                bottom_xvals = [x*-aSign for x in bottom_xvals]
                bottom_xvals.sort()
                bottom_yvals = [y] * len(bottom_xvals)
                bottom_zvals = [z] * len(bottom_xvals)
                Height = ExBankHeight + exZadj + (y*slopeDiff)
                TopWidth = ((Height*math.tan((90-EqBankAngle)*math.pi/180))*2)+ExBottomWidth
                bank_xvals = list(range(math.ceil(ExBottomWidth/2),math.ceil(TopWidth/2)))
                bank_xvals.extend([TopWidth/2])
                bank_xvals = [x*-aSign for x in bank_xvals]
                bank_xvals.sort()
                # the corners to interpolate between need to be re-ordered based on aSign
                bank_interpX = [ExBottomWidth/2*-aSign,TopWidth/2*-aSign] if aSign<0 else [TopWidth/2*-aSign,ExBottomWidth/2*-aSign]
                bank_interpY = [(y*EqBedSlope),Height+(y*EqBedSlope)] if aSign<0 else [Height+(y*EqBedSlope),(y*EqBedSlope)]
                bank_zvals = list(np.interp(bank_xvals,
                                    bank_interpX,bank_interpY))
                bank_yvals = [y] * len(bank_xvals)
                
                # and rotate this every X degrees as needed.
                aRngInc  = int(aDiffInc*aSign)
                aRangEnd = int(math.ceil(abs(aDiff))*aSign)
                for aExtra in range(aRngInc,aRangEnd,aRngInc):
                    # make df for each extra XS we want to rotate
                    dfExtra = pd.DataFrame({'X': bottom_xvals+bank_xvals,'Y': bottom_yvals+bank_yvals})
                    rotExtra = rotate(dfExtra, origin=(0,0), degrees=aExtra)
                    Xlist.extend(list(rotExtra[0]))
                    Ylist.extend(list(rotExtra[1]))
                    Zlist.extend(bottom_zvals+bank_zvals)
            
            # make dataframe from X and Y lists
            dfxy = pd.DataFrame({'X': Xlist,'Y': Ylist})
            
            # Shift everything to our starting point XY (no longer [0,0] but state plane coordinates)
            dfxy['X'] = dfxy['X'] + stX
            dfxy['Y'] = dfxy['Y'] + stY

            # Rotate the point cloud around starting point, based on the bearing of the line segment
            Rot = rotate(dfxy, origin=(stX,stY), degrees=-angle)

            # Make it a dataframe, add Z back in, and shift elevation to match starting Z
            Rotxyz = pd.DataFrame(Rot)
            Rotxyz = Rotxyz.rename(columns={0:'X', 1:'Y'})
            Rotxyz['Z'] = Zlist
            Rotxyz['Z'] = Rotxyz['Z'] + stZ
            
            # return point cloud array
            arr = pd.DataFrame(Rotxyz,columns=['X','Y','Z']).to_records(index=False)
            return arr
            
        ## End functions ##
        
        # determine the number of input features and segments. Was going to use this for a progressor bar but not anymore.
        numSegList = []
        features = [f[0] for f in arcpy.da.SearchCursor('line_lyr',"SHAPE@")]
        for f in features:
            numSegList.append(f.pointCount-1)
        messages.addMessage('Processing {} input line features containing {} total segments'.format(len(oids),sum(numSegList)))
        
        # If flip check is on, check line direction and reverse if necessary
        if flip_check:
            messages.addMessage('Checking for downhill lines...')
            # get endpoints of lines
            tempPt = arcpy.CreateUniqueName("ptDir", arcpy.env.scratchGDB)
            arcpy.management.GeneratePointsAlongLines('line_lyr', tempPt, "PERCENTAGE", None, 100, "END_POINTS")
            
            # assign Z value
            arcpy.sa.AddSurfaceInformation(tempPt, in_raster_ref, "Z", "BILINEAR", None, 1, 0, '')

            # Loop through ORIG_FID of points, and compare first point to last point
            orig_oids = set(row[0] for row in arcpy.da.SearchCursor(tempPt, "ORIG_FID"))
            flip_oids = []
            for orig in orig_oids:
                # get z values for that OID
                lineZs = [row[0] for row in arcpy.da.SearchCursor(tempPt,'Z',where_clause = 'ORIG_FID = '+str(orig))]
                # if start is higher than end, need to flip direction
                if lineZs[0] > lineZs[1]:
                    flip_oids.append(orig)
                    messages.addMessage('--OID {} is downhill and will be flipped'.format(orig))
            # If there are features to flip, select them and flip them
            if len(flip_oids) > 0:
                messages.addMessage('Reversing direction for {} feature(s)...'.format(len(flip_oids)))
                arcpy.management.MakeFeatureLayer('line_lyr', 'flip_lyr', 
                                            where_clause=arcpy.Describe(in_line).OIDFieldName+" IN ({:s})".format(','.join(f"{x}" for x in flip_oids)))
                arcpy.edit.FlipLine('flip_lyr')
            else:
                messages.addMessage('No downhill lines found.')
            arcpy.management.Delete(tempPt)
        
        
        ## Process the input line file
        # Make sure we have minimum elevation for each line.
        arcpy.sa.AddSurfaceInformation(in_line, in_raster_ref, "Z_MIN;Z_MAX", "BILINEAR", None, 1, 0, '')
                    
        # get lines as numpy array
        arr_fields = ['OID@','SHAPE@XY','Z_Min','Z_Max','bedType', 'drainageAcres', 'exBedSlope', 'manningsN', 'flowRate', 'initialValDepth', 
                        'exBottomWidth', 'exBankSlope', 'shieldsParam', 'critBedParam', 'chanFormDischperUnitWidth', 'meanGrainSize', 
                        'designDisch', 'medGrainSize','eqBedSlope','eqBankSlope', 'exBankHeight']
        line_arr = arcpy.da.FeatureClassToNumPyArray('line_lyr', arr_fields, explode_to_points=True)
        
        # get dictionary of DEM Z_Min and OID
        d_z_demInit = {row[0]: row[1] for row in arcpy.da.SearchCursor('line_lyr',[arcpy.Describe('line_lyr').OIDFieldName,'Z_MIN'])}
        
        # also initialize dictionaries of final Eq and Ex bed elevations. So we have something to pull from for non-dependent reaches
        # these can just be copies of d_z_demInit to begin with
        d_z_eq  = d_z_demInit.copy()
        d_z_ex  = d_z_demInit.copy()
        
        # also need a dict to lookup fromOID from OID
        if start_Z == 'Continuous (link using fromOID)':
            # if fromOID is empty use OID
            d_oid  = {row[0]: row[0] if row[1] is None else row[1] for row in arcpy.da.SearchCursor('line_lyr',[arcpy.Describe('line_lyr').OIDFieldName,'fromOID'])}
        
        # prepare empty lists of poly bounding boxes and rasters (for each segment to be compiled at the end), and a list of items to delete
        del_list  = []
        poly_list = []
        ras_list  = []
        
        # If using method 'Link using fromOID', we need to process features in a particular order.
        if start_Z == 'Continuous (link using fromOID)':
            # get tuple of (OID,fromOID)
            listDepend = [(row[0],row[1]) for row in arcpy.da.SearchCursor('line_lyr',[arcpy.Describe('line_lyr').OIDFieldName,'fromOID'])]
            # if OID = fromOID, set fromOID to None
            listDepend = [(item,None if dep==item else dep) for item,dep in listDepend]
            
            # Use topological sorting to order the OIDs
            # a simplified Kahn's algo (https://en.wikipedia.org/wiki/Topological_sorting)
            oid_order = []
            s = [item for item,dep in listDepend if dep is None]
            while len(s) > 0:
                n = s.pop(0)
                oid_order.append(n)
                for m in [item for item,dep in listDepend if dep is n]:
                    s.append(m)
        else:
            oid_order = set(line_arr['OID@']) 
        
        # loop through each OID
        for oid in oid_order:
            # extract array for just that feature
            feat_arr = line_arr[line_arr['OID@']==oid]
            numSeg = len(feat_arr)-1
            
            # Calculate exBedSlope from DEM if not defined (assumes everything is in feet)
            if np.isnan(feat_arr['exBedSlope'][0]):
                # slope is change in Z over length of this feature
                feat_arr['exBedSlope'][0] = (feat_arr['Z_Max'][0]-feat_arr['Z_Min'][0]) / [r[0] for r in arcpy.da.SearchCursor('line_lyr','SHAPE@LENGTH',where_clause = arcpy.Describe('line_lyr').OIDFieldName+"="+str(oid))][0]
                # update the input FC
                with arcpy.da.UpdateCursor('line_lyr', 'exBedSlope',where_clause = arcpy.Describe('line_lyr').OIDFieldName+"="+str(oid)) as cursor:
                    for row in cursor:
                        row[0] = feat_arr['exBedSlope'][0]
                        cursor.updateRow(row)
            
            # calculate EqBedSlope (if value is missing, or box is checked), update feat_arr
            if calc_eq or np.isnan(feat_arr['eqBedSlope'][0]):
                if feat_arr['bedType'][0] == 'Cohesive':
                    slopeCalc = cohesive_EqBedSlope(feat_arr['drainageAcres'][0])
                    feat_arr['eqBedSlope'][0] = slopeCalc
                elif feat_arr['bedType'][0] == 'Sand and Gravel': 
                    depth = normal_depth(feat_arr['exBedSlope'][0], feat_arr['manningsN'][0], feat_arr['flowRate'][0], feat_arr['initialValDepth'][0], feat_arr['exBottomWidth'][0], feat_arr['exBankSlope'][0])
                    slopeCalc = sand_gravel_EqBedSlope(depth)
                    feat_arr['eqBedSlope'][0] = slopeCalc
                elif feat_arr['bedType'][0] == 'Coarser Than Sand':
                    slopeCalc = coarserThanSand_EqBedSlope(feat_arr['shieldsParam'][0], feat_arr['critBedParam'][0], feat_arr['chanFormDischperUnitWidth'][0], feat_arr['manningsN'][0], feat_arr['meanGrainSize'][0], feat_arr['designDisch'][0], feat_arr['medGrainSize'][0])
                    feat_arr['eqBedSlope'][0] = slopeCalc
                else: # "Custom" bedType for manual definition of eq slope
                    slopeCalc = feat_arr['eqBedSlope'][0]
                
                # Update eqBedSlope on the input line data as needed
                with arcpy.da.UpdateCursor('line_lyr', 'eqBedSlope',where_clause = arcpy.Describe('line_lyr').OIDFieldName+"="+str(oid)) as cursor:
                    for row in cursor:
                        row[0] = slopeCalc
                        cursor.updateRow(row)
            
            # extract parameters for that feature
            ExBedSlope    = feat_arr['exBedSlope'][0]
            EqBedSlope    = feat_arr['eqBedSlope'][0]
            slopeDiff     = ExBedSlope-EqBedSlope
            EqBankSlope   = feat_arr['eqBankSlope'][0]
            ExBottomWidth = feat_arr['exBottomWidth'][0]
            ExBankHeight  = feat_arr['exBankHeight'][0]
            
            
            # get elevations depending on method
            if start_Z == 'Continuous (link using fromOID)':
                # get initial elevation for feature
                stZ = d_z_eq[d_oid[oid]]
                # and track bed elevation of existing channel
                exZ = d_z_ex[d_oid[oid]]
            else:
                stZ = d_z_demInit[oid]
                exZ = d_z_demInit[oid]
            
            exZadj = exZ - stZ
            
            arcpy.AddMessage('--OID {}: Contains {} segments. Starting bed elevation {:.3f} feet'.format(oid,len(feat_arr)-1,stZ))
            arcpy.AddMessage('---Equilibrium Bed Slope: {:.5f}'.format(EqBedSlope))
            
            ## loop through each segment within the feature
            # because we're looking for sharp corners in the channel with lastAngle, prepopulate here for use in the first segment
            lastAngle = getBearing(feat_arr['SHAPE@XY'][0][0],feat_arr['SHAPE@XY'][0][1],feat_arr['SHAPE@XY'][1][0],feat_arr['SHAPE@XY'][1][1])

            for ptNum in range(numSeg):
                # get the starting/ending coordinates of the segment
                stX = feat_arr['SHAPE@XY'][ptNum][0]
                stY = feat_arr['SHAPE@XY'][ptNum][1]
                enX = feat_arr['SHAPE@XY'][ptNum+1][0]
                enY = feat_arr['SHAPE@XY'][ptNum+1][1]
                angle = getBearing(stX,stY,enX,enY)
                # get the point cloud for that segment
                seg_arr = ReturnPointCloudArray(EqBedSlope,EqBankSlope,ExBottomWidth,ExBankHeight,stX,stY,enX,enY,stZ,exZadj,slopeDiff,angle,lastAngle)
                lastAngle = angle
                # get distance of that segment to track eq and ex Zs
                segDist = math.sqrt((enX - stX)**2 + (enY - stY)**2)
                enZ     = stZ + (segDist*EqBedSlope)
                exZ     = exZ + (segDist*ExBedSlope)
                # next segment stZ is this segment's enZ
                stZ = enZ
                # and the difference is our Zadj for the next segment
                exZadj = exZ - enZ
                
                arcpy.AddMessage('----Segment {}: Distance: {:06.2f}, Eq rise: {:.3f}, Ex rise: {:.3f}, Angle: {:.0f}'.format(ptNum+1,segDist,(segDist*EqBedSlope),(segDist*ExBedSlope),angle))
                
                # compile an array of points for the feature from the segments
                if ptNum == 0:
                    point_arr = seg_arr
                else:
                    point_arr = np.concatenate((point_arr, seg_arr))
                
                # each segment's point cloud needs it's own bounding geometry, which we can dissolve to get the channel polygon          
                tempPt = arcpy.CreateUniqueName("ptSeg", arcpy.env.scratchGDB)
                arcpy.da.NumPyArrayToFeatureClass(seg_arr, tempPt, ['X','Y','Z'], arcpy.SpatialReference(2893))
                tempPoly = arcpy.CreateUniqueName("polySeg", arcpy.env.scratchGDB)
                arcpy.management.MinimumBoundingGeometry(tempPt, tempPoly, 'CONVEX_HULL')
                poly_list.append(tempPoly)
                del_list.append(tempPoly)
                
                # Previous versions did the IDW for each feature, but we got weird noisy behavior for sharp bends at the top of steep channels
                # this is for less efficient, but gives a better product. Create IDW for each segment, and mosaic at the end
                with arcpy.EnvManager(snapRaster=in_raster_ref, mask=tempPoly):
                    tempRas = arcpy.sa.Idw(tempPt, "Shape.Z", cell_size=1)
                ras_list.append(tempRas)
                del_list.append(tempRas)
                arcpy.management.Delete(tempPt)
                
            #Update the dictionaries for the feature
            d_z_eq[oid] = enZ
            d_z_ex[oid] = exZ
            
            #print the final bed elevation for the feature
            messages.addMessage('----Final bed elevation: {:.3f} feet:'.format(enZ))
            
        # end OID loop     
        
        # merge all the segment polygons together
        messages.addMessage('Merging segment polygons...')
        tempPolyMerge = arcpy.CreateUniqueName("polyM", arcpy.env.scratchGDB)
        arcpy.management.Merge(poly_list, tempPolyMerge)
        del_list.append(tempPolyMerge)
        # and dissolve to remove overlaps
        tempPolyDiss = arcpy.CreateUniqueName("polyD", arcpy.env.scratchGDB)
        arcpy.management.Dissolve(tempPolyMerge, tempPolyDiss)
        del_list.append(tempPolyDiss)
        # We see small slivers missing from these polygons. A 5 foot buffer, reversed, will fill them
        tempPolyBuff = arcpy.CreateUniqueName("poly", arcpy.env.scratchGDB)
        arcpy.analysis.Buffer(tempPolyDiss, tempPolyBuff, "5 Feet", "FULL", "ROUND", "NONE", None, "PLANAR")
        del_list.append(tempPolyBuff)
        # output polygon
        arcpy.analysis.Buffer(tempPolyBuff, out_poly, "-5 Feet", "FULL", "ROUND", "NONE", None, "PLANAR")
        arcpy.management.DeleteField(out_poly, 'BUFF_DIST')
        arcpy.management.DeleteField(out_poly, 'ORIG_FID')
        
        # mosaic the segment rasters from all features, and mask to out_poly 
        messages.addMessage('Combining segment rasters...')
        tempMosaic = arcpy.CreateUniqueName("mosaic", arcpy.env.scratchGDB)
        del_list.append(tempMosaic)
        with arcpy.EnvManager(snapRaster=in_raster_ref, mask=out_poly):
            arcpy.management.MosaicToNewRaster(ras_list, os.path.dirname(tempMosaic), os.path.basename(tempMosaic), 
                                            coordinate_system_for_the_raster = arcpy.SpatialReference(2893),
                                            pixel_type='32_BIT_FLOAT', cellsize=1, number_of_bands=1, 
                                            mosaic_method='MINIMUM')
        tempExtract = arcpy.sa.ExtractByMask(tempMosaic, out_poly); 
        del_list.append(tempExtract)
        
        ## fill in any missing data slivers or gaps
        # To achieve this, convert raster to points, then IDW interp again for the final polygon
        messages.addMessage('Preparing final raster...')
        tempPt = arcpy.CreateUniqueName("ptFromRas", arcpy.env.scratchGDB)
        del_list.append(tempPt)
        arcpy.conversion.RasterToPoint(tempExtract, tempPt)
        tempIDW     = arcpy.CreateUniqueName("tempIDW", arcpy.env.scratchGDB)
        del_list.append(tempIDW)
        with arcpy.EnvManager(snapRaster=in_raster_ref, mask=out_poly):
            # for some reason, SA idw doesn't let you set output raster name, but 3d idw does.
            arcpy.ddd.Idw(tempPt, "grid_code", tempIDW, cell_size=1)
            
        ## Zonal stats - erosion credit over 30 years (50% efficiency)
        # Project / Resample input DEM to 1 foot grid, snapped and masked to our final raster output
        messages.addMessage('Computing difference raster...')
        tempDEM     = arcpy.CreateUniqueName("tempDEM", arcpy.env.scratchGDB)
        tempMinus   = arcpy.CreateUniqueName("tempMinus", arcpy.env.scratchGDB)
        del_list.append(tempDEM)
        del_list.append(tempMinus)
        with arcpy.EnvManager(snapRaster = tempIDW, cellSize = 1, extent = out_poly, mask = out_poly):
            arcpy.management.ProjectRaster(in_raster_ref, tempDEM, arcpy.SpatialReference(2893), 'BILINEAR', "1")
            #Get minimum value of the two rasters (to prevent modeled elevation being above existing elevation)
            min_raster = Min([tempIDW,tempDEM])
            min_raster.save(out_ras)
            # Subtract eq raster from this DEM
            arcpy.ddd.Minus(tempDEM, out_ras, tempMinus)
        
        # Zonal stats
        messages.addMessage('Adding Zonal Statistics to output polygon...')
        tempTbl = arcpy.CreateUniqueName("tempTbl", arcpy.env.scratchGDB)
        del_list.append(tempTbl)
        ZonalStatisticsAsTable(out_poly, arcpy.Describe(out_poly).OIDFieldName, tempMinus, tempTbl,"DATA","SUM")
        arcpy.management.JoinField(out_poly, arcpy.Describe(out_poly).OIDFieldName, tempTbl, arcpy.Describe(out_poly).OIDFieldName, ["SUM"])
        # Add and populate field in output polygon
        arcpy.management.AddField(out_poly, "Vol_yr", "DOUBLE")
        arcpy.management.CalculateField(out_poly, "Vol_yr", "!SUM! / 60","PYTHON3")
        arcpy.management.DeleteField(out_poly,"ORIG_FID")
        arcpy.management.DeleteField(out_poly,"SUM")
        
        # delete intermediate files
        messages.addMessage('Deleting intermediate files...')
        arcpy.management.Delete(del_list)
        #arcpy.management.Delete('line_lyr')
        messages.addMessage('Script Complete!')
        
        return



class setRasStretch(object):

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Set Raster Symbology"
        self.description = "Quickly set multiple rasters to custom min/max stretch symbology"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        
        rasList = arcpy.Parameter(
            displayName="Input Raster(s)",
            name="rasList",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input",
            multiValue = True)
            
        minVal = arcpy.Parameter(
            displayName="Minimum stretch value",
            name="minVal",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        minVal.value = 0
        
        maxVal = arcpy.Parameter(
            displayName="Maximum stretch value",
            name="maxVal",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        maxVal.value = 200
        
        cRamp = arcpy.Parameter(
            displayName="Define colorramp",
            name="cRamp",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        
        p = arcpy.mp.ArcGISProject('CURRENT')
        # Get elevation color ramps, sort them
        cList = [c.name for c in p.listColorRamps('Elev*')]
        cSort = [i for (v, i) in sorted((v, i) for (i, v) in enumerate([int(x.split("#")[-1]) for x in cList]))]
        cList = [cList[i] for i in cSort]
        cRamp.filter.list = cList
        cRamp.value = 'Elevation #1'

        params = [rasList,minVal,maxVal,cRamp]
        
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        
        #Get the input data
        rasList =  parameters[0].values
        minVal  =  parameters[1].Value
        maxVal  =  parameters[2].Value
        cRamp   =  parameters[3].valueAsText
        
        #get current project and map view
        p = arcpy.mp.ArcGISProject('CURRENT')
        m = p.activeMap

        for ras in rasList:
            messages.addMessage('Updating symbology for raster: {}'.format(ras.name))
            #get raster layer
            lyr  = m.listLayers(ras.name)[0]

            #set symbology (color ramp, labels, stretch type)
            sym = lyr.symbology
            sym.updateColorizer('RasterStretchColorizer')
            sym.colorizer.stretchType = 'MinimumMaximum'
            cr = p.listColorRamps(cRamp)[0]
            sym.colorizer.colorRamp = cr
            sym.colorizer.minLabel = "Min: " + str(minVal)
            sym.colorizer.maxLabel = "Max: " + str(maxVal)
            lyr.symbology = sym

            #use CIM to set custom statistics
            cim_lyr = lyr.getDefinition('V2')
            #cim_lyr.colorizer.statsType = 'GlobalStats'
            cim_lyr.colorizer.useCustomStretchMinMax = True
            cim_lyr.colorizer.customStretchMin = minVal
            cim_lyr.colorizer.customStretchMax = maxVal
            #cim_lyr.colorizer.stretchStats.min = minVal
            #cim_lyr.colorizer.stretchStats.max = maxVal
            lyr.setDefinition(cim_lyr)
        
        return
