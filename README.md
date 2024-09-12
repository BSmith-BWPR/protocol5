# Repository
Python toolbox for Protocol 5 Crediting for outfall gully restoration projects

Included here is the Python Toolbox and the file geodatabase (zipped) containing an empty line feature class template for use with the tool.

Also included are PDFs of the Expert Panel and MDOT documents for Protocol 5.

# Protocol 5 Tool
## Background

The Expert Panel document Recommendations for Crediting Outfall and Gully Stabilization Projects in the Chesapeake Bay Watershed defines the crediting mechanism of Protocol 5 (Alternative Prevented Sediment for Outfalls).  This approach accounts for sediment loss through vertical incision at outfalls and gullies.  Credit is estimated by comparing existing channel conditions to modeled equilibrium channel conditions. For more information on the Protocol, please review the [Expert Panel guidance](https://www.chesapeakebay.net/channel_files/37043/approval_draft_outfall_restoration_memo_070119.pdf). 

This python toolbox produces GIS raster and polygon data for the modeled equilibrium channel condition.  Total sediment volume is converted to an annual sediment volume credit, assuming 50% efficiency and a 30 year time window as outlined in the Expert Panel document.


## Requirements

- ArcGIS Pro 2.x
-	Spatial Analyst and 3D Analyst extensions 
-	Python


## Use

The tool is designed to be run for a single restoration project, but a project can contain multiple line features with different attributes.  Processing time is dependent on the number of vertices in the input line feature, so using the fewest vertices necessary to accurately capture channel geometry will shorten computation time.

Line segments should be drawn uphill – beginning at the downstream grade control and ending at the upstream extent.  The tool will reverse any downhill lines if the box “Check for / enforce uphill line direction” is selected.  This box is unselected by default to reduce processing time.

DEM vertical units must be in Feet. Map units must be in Feet.  The preferred spatial reference is WKID 2893 (NAD 1983 HARN StatePlane Maryland FIPS 1900 US Feet). Output rasters and feature classes will use this spatial reference.

Slope inputs should be in V:H format (where lower values correspond to shallower slopes).


## Geodatabase with Template Line Feature Class

The ‘P5LineTemplate.gdb’ geodatabase contains a line feature class with all the required fields.  The following fields must be populated for each feature:
- exBottomWidth:	Existing channel bottom width (feet)
- exBankHeight:	Existing channel bank height (feet)
- exBedSlope:	Existing channel bed slope (feet)
- bedType:	Bed Material.  (Custom, Cohesive, Sand and Gravel, Coarser than Sand)
- eqBankSlope:	Equilibrium bank slope (V:H), Default 0.568

## Estimating Equilibrium Slope Based On Bed Material
Equilibrium slope is estimated from channel properties.  Depending on the bed material type, different parameters are required to estimate equilibrium slope.  The four options for bed type, and associated required fields, are shown below.
1. Custom
    - eqBedSlope:	Equilibrium channel bed slope, entered directly by the user

2. Cohesive
    - drainageAcres:	The contributing drainage area to the reach, in acres

3. Sand and Gravel
    - exBankSlope:	Existing bank slope / channel side slope (V:H)
    - manningsN:	Manning’s roughness coefficient n, default 0.025
    - flowRate:	10-year flow rate (CFS)
    - initialValDepth:	Initial estimate of normal depth (feet), default 0.2

4. Coarser Than Sand
    - manningsN:	Manning’s roughness coefficient, n, default 0.025
    - shieldsParam:	Shield’s parameter, θc
    - critBedParam:	Critical bed material size, Dc (feet)
    - chanFormDischperUnitWidth:	Channel forming discharge per unit width, q (feet2/second)
    - designDisch:	Design discharge, Qd (CFS)
    - meanGrainSize:	Mean grain size Dm (mm)
    - medGrainSize:	Median grain size D50 (mm)

See the MDOT SHA document Alternative Headwater Channel and Outfall Crediting Protocol for additional information.


## Reach Connectivity – Continuous Equilibrium Elevation
A project may contain reaches with different channel attributes that must be represented by multiple connected line features.  The tool parameter Reach Connectivity Method has two options:
1. Segmented – Reaches are processed in distinct segments, using the DEM for initial elevation
2. Continuous – Reaches are processed as a connected channel, where the initial elevation is determined by the final elevation of the downstream connected reach

![image](https://user-images.githubusercontent.com/103444219/173117388-b73a76e2-83df-45bf-a954-0124318f59a0.png)


If using the continuous method, the fromOID field in the template line feature class must be populated.  In the example above, the fromOID value of Reach 3 would be the ObjectID of Reach 2, and the fromOID value of Reach 2 would be the ObjectID of Reach 1.  A feature does not need a value in the fromOID field if its initial elevation is determined by the DEM and not another feature (in the example above, fromOID for Reach 1 would be blank). 

## Output
The tool generates raster and polygon feature class output.  The default output workspace is the Scratch Geodatabase environment variable.

The raster output is in feet, with a cell size of 1 foot.

The polygon feature class is the extent of the equilibrium channel raster. It contains a field Vol_yr containing the annual credit volume in cubic feet, assuming a 30 year period and 50% efficiency.

# Set Raster Symbology Tool
The Python toolbox also includes the Set Raster Symbology Tool.  This tool applies a uniform symbology to multiple raster layers with a single operation:

- Primary Symbology: 	Stretch
- Stretch Type: 		Minimum Maximum
- Statistics: 		Custom
- Min / Max values: 	User Defined
- Color Ramp: 		User Defined

Note that due to the limitations of the python CIM functions, the changes made to the ‘Statistics’ tab in the Symbology Window will not appear until closing and reopening the project, even though the raster symbology itself will immediately update.

# Disclaimer
By acceptance of the GIS material, you agree as follows: 
The GIS material (the “material”) is made available by Anne Arundel County, Maryland (the “County”) as a public service. The material is for reference purposes only, and the County makes no representatives, warranties, or guarantees of the accuracy of the material. THE COUNTY MAKES NO AND DISCLAIMS ALL EXPRESS AND IMPLIED WARRANTIES RELATING TO THE MATERIAL, INCLUDING WARRANTIES OF MERCHANTABILITY, INTEGRATION, TITLE, AND FITNESS FOR A PARTICULAR PURPOSE. You release the County, its agents, servants, and employees, from any and all liability related to the material or any of it, including its accuracy, availability, use, and misuse. In no event shall the County be liable for any direct, indirect, incidental, consequential, or other damages, including savings, profits, fees, costs, loss of data, or business interruption, related in any way to the material or any of it, including its accuracy, availability, use, and misuse. 
The material is in the public domain and may be copied without permission. Citation to the source is requested. 
Any errors or omissions in the material should be reported to the Anne Arundel County Bureau of Watershed Protection and Restoration, TMDL Support Group. 
