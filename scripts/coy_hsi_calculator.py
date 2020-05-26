import arcpy
from arcpy.sa import *


# Define Functions-------------------------------------------------------------
def calc_edge_area(outer_raster,inner_raster,radius=200):
    """
    Calculates an "Edge" raster; an area of interest (e.g. forest) 
    that edges a cover type of interest (e.g. wetlands, streets, etc). 
    
    Args:
        outer_raster: The habitat of interest that will make up the edge 
                      environment. This raster should be boolean, with
                      cells that are of the cover type of interest
                      (e.g. forest) having a score of 1 and all other cells
                      scored 0
        inner_raster: The inner habitat that the outer raster edges. Should
                      be a boolean raster where cells of interest are scored 1
                      and all others 0
        radius: Controls the width of the edge raster returned, by default 
                200 map units will be used.
                
    Returns:
        A raster of the edge area defined per the provided arguments

    """
    
    # Run focal stats tool to find cells that "edge" the input cover_type
    # output raster will look similar to the original, but cells that edge the
    # cover will now have a score of 1 instead of 0                                  
    inner_edge = FocalStatistics(inner_raster,                 
                                        NbrRectangle(3, 3, "CELL"), 
                                        "MAXIMUM", 
                                        True) 
    
    # find cells in the outer edge raster that are adjacent the 
    # inner edge
    outer_cells_adj_inner = ((outer_raster == 1)&(inner_edge > 0))
    
    # Run focal stats to find cells within x map unit radius
    # of the boundary betweent he two habitat rasters
    nbr_circle = NbrCircle(radius,"MAP")
    within_x_of_edge = FocalStatistics(outer_cells_adj_inner, 
                                       nbr_circle, 
                                       statistics_type = "MAXIMUM", 
                                       ignore_nodata = True)
    
    
    # Create new raster where cells of the outer habitat are of the cover type
    # of interest (e.g. forest) and within x radius of the habitat boundary
    edge_raster = ((outer_raster == 1)&(within_x_of_edge > 0))
    
    return edge_raster


def reclassify_wrapper(raster,remap_obj,field = "Value"):
    """
    Wrapper around the reclassify tool
    
    Args:
        raster: The raster to be reclassified
        remap_obj: a RemapValue or RemapRange object
        field: Field to act on, "Value" by default.
        
    Returns:
        A reclassified raster layer
  
    """
    
    reclass = Reclassify(raster,field,remap_obj)
    
    return reclass

    
def rescale_wrapper(dist_layer,function,midpoint="#",
                    spread="#",scale_start=1,scale_end=100):
    """
    Wrapper around the rescale tool that limits the transformation
    functions used to "small" and "large". Only dist_layer and function
    are required; midpoint will default to the mean of the dist_layer, and
    spread to a value of 5.
    
    Args:
        dist_layer: The a raster layer that represents a disturbance factor.
        function: Transformation function to use (must be "small" or "large").
        midpoint: Values greater than midpoint decrease in preference("small").
                  Values greater than midpoint increase in preference("large").
        spread: How quickly values increase and decrease as they move from the
                midpoint.
        scale_start: The start of the evaluation scale (e.g. 1).
        scale_end: The end of the evaluation scale (e.g. 100).
    
    Returns:
        a rescaled raster layer
        
    Raises:
        Exception: Raised if non-valid function parameter is supplied
    
    """
        
    if function.lower() == "small":    
        transformation_function = TfSmall(midpoint, spread)
    elif function.lower() == "large":
        transformation_function = TfLarge(midpoint, spread)
    else:
        raise Exception("Must use either the 'small' or 'large' function, "
                        "check documentation and try again")

    reclass = RescaleByFunction(dist_layer,transformation_function, 
                                scale_start,scale_end)

    return reclass

    
def calc_building_factor(buildings,cellsize,radius,midpoint="#",
                         spread="#",scale_start=1,scale_end=100):
    """
    Calculates the building factor for the Eastern Coyote Habitat Suitability
    Model. Building points per the specified radius are calculated to gauge
    building density and the a rescale is performed.
    
    Args:
        buildings: Point or polygon features that represent buildings.
        radius: Radius used to calculate building density per unit area
        midpoint: Values greater than midpoint decrease in preference.
        spread: How quickly values increase and decrease as they move from the
                midpoint.
        scale_start: The start of the evaluation scale (e.g. 1).
        scale_end: The end of the evaluation scale (e.g. 100).
    
    Returns:
        The building factor for the model per the specified arguments
        
    Raises:
        Exception: Raised if a feature class of the wrong shapetype is provided
    
    """
    shape_type = arcpy.Describe(buildings).shapeType                                             
    if shape_type == 'Polygon':
        building_pts = arcpy.FeatureToPoint_management(buildings, 
                                                       r"in_memory\buildPts", 
                                                       "INSIDE")
    elif shape_type == 'Point':
        building_pts = buildings
        
    else:
        raise Exception("buildings must be polygon or point!")
        
        
    pt_density = PointDensity(building_pts, "NONE",cellsize, 
                              NbrCircle(radius,"MAP"))
    
    factor = rescale_wrapper(pt_density,"small",midpoint,
                             spread,scale_start,scale_end)
    
    arcpy.Delete_management(r"in_memory\buildPts")
                                
    return factor


def calc_street_factor(streets,cellsize,midpoint="#",
                       spread="#",scale_start=1,scale_end=100):
    """
    Calculates the street factor for the Eastern Coyote Habitat Suitability
    Model.Euclidean distance is calculated for cells in the raster, further
    a cell is from a street the higher it will be scored
    
    Args:
        streets: line feature representing streets/roads.
        cellsize: Cell size to use for euclidean distance calculation
        midpoint: Values greater than midpoint increase in preference.
        spread: How quickly values increase and decrease as they move from the
                midpoint.
        scale_start: The start of the evaluation scale (e.g. 1).
        scale_end: The end of the evaluation scale (e.g. 100).
    
    Returns:
        The building factor for the model per the specified arguments
        

    """

    dt_streets = EucDistance(streets,cell_size=cellsize)
    
    factor = rescale_wrapper(dt_streets,"large",midpoint,
                             spread,scale_start,scale_end)
    
    return factor


def calculate_hsi(lc_ras,forest_ras,edge_width,wetland_ras,ag_ras,lc_remap,
                  canopy_ras,canopy_remap,buildings,streets,building_radius,
                  building_mp='#',building_spread='#',street_mp='#',
                  street_spread='#',scale_start=1,scale_end=100):
    """
    Calculates a Habitat Suitability Index for the Eastern coyote
    
    Args:
        lc_ras: A land-cover raster
        forest_ras: forested areas derived from the provided lc_ras
        edge_width: Desired width of the edge areas calculated for the model
          default is 200 map units
        wetland_ras: same as forest_Ras but wetlands
        ag_ras: same as above but agriculture/crop/pasture lands, etc.
        lc_remap: A remap value or remap range object
        canopy_ras: A raster representing canopy coverage for the aoi
        canopy_remap: A remap value or remap range object
        buildings: A feature class representing either building points or
                   polygons
        streets: A feature class representing street lines
        building_mp: Values greater than midpoint decrease in preference.
        building_spread: How quickly values increase and decrease as they 
                         move from the midpoint.
        street_mp: same as building_mp
        street_spread: same as building_spread
        scale_start: The start of the evaluation scale (e.g. 1)
        scale_end: The end of the evaluation scale (e.g. 100)
        
    Returns:
        A raster layer representing HSI scores for the Eastern
    
        
    """
    
    edge_remap = RemapValue([[1,scale_end],[0,scale_start]])
        
    for_edge_wet_factor = reclassify_wrapper(calc_edge_area(forest_ras,
                                                            wetland_ras,
                                                            edge_width),
                                             edge_remap)
    
    for_edge_ag_factor = reclassify_wrapper(calc_edge_area(forest_ras,
                                                           ag_ras,
                                                           edge_width),
                                            edge_remap)
    
    building_factor = calc_building_factor(buildings,lc_ras,building_radius,
                                           building_mp,building_spread,
                                           scale_start,scale_end)
    
    street_factor = calc_street_factor(streets,lc_ras,
                                       street_mp,street_spread,
                                       scale_start,scale_end)
    
    lc_factor = reclassify_wrapper(lc_ras,lc_remap)
    
    canopy_factor = reclassify_wrapper(canopy_ras,canopy_remap)
    
    factors = [for_edge_wet_factor,for_edge_ag_factor,building_factor,
               street_factor,lc_factor,canopy_factor]
    
    # Ensures that resulting raster edges remain intact
    # converts no data values to 0
    factors = [Con(IsNull(factor),0,factor) for factor in factors]

    hsi = sum(factors)/6
    
    final_result = [hsi] + factors
    
    return final_result
    
    
if __name__=='__main__':
     
    aoi = arcpy.GetParameter(0)
    wks = arcpy.GetParameterAsText(1)
    save_factors = arcpy.GetParameter(2)
    lc_ras = arcpy.Raster(arcpy.GetParameterAsText(3))
    lc_remap = arcpy.GetParameterAsText(4)
    forest_ras = arcpy.Raster(arcpy.GetParameterAsText(5))
    edge_width = int(arcpy.GetParameterAsText(6))
    wetland_ras = arcpy.Raster(arcpy.GetParameterAsText(7))
    ag_ras = arcpy.Raster(arcpy.GetParameterAsText(8))
    canopy_ras = arcpy.Raster(arcpy.GetParameterAsText(9))
    canopy_remap = arcpy.GetParameterAsText(10)
    buildings = arcpy.GetParameterAsText(11)
    streets = arcpy.GetParameterAsText(12)
    building_radius = int(arcpy.GetParameterAsText(13))
    # deal with optional parameters
    try:
        building_mp = int(arcpy.GetParameterAsText(14))
    except:
        building_mp = '#'
    try:
        building_spread = int(arcpy.GetParameterAsText(15))
    except:
        building_spread = '#'
    try:
        street_mp = int(arcpy.GetParameterAsText(16))
    except:
        street_mp = '#'
    try:
        street_spread = int(arcpy.GetParameterAsText(17))
    except:
        street_spread = '#'
        
    scale_start = int(arcpy.GetParameterAsText(18))
    scale_end = int(arcpy.GetParameterAsText(19))

    arcpy.env.overwriteOutput = True
    arcpy.env.extent = (aoi)

    
    result = calculate_hsi(lc_ras,forest_ras,edge_width,wetland_ras,ag_ras,
                           lc_remap,canopy_ras,canopy_remap,buildings,
                           streets,building_radius,building_mp,building_spread,
                           street_mp,street_spread,scale_start,scale_end)
    
    if save_factors:
        result[1].save(wks+r"\for_edge_wet_factor")
        result[2].save(wks+r"\for_edge_ag_factor")
        result[3].save(wks+r"\building_factor")
        result[4].save(wks+r"\lc_factor")
        result[5].save(wks+r"\street_factor")
        result[6].save(wks+r"\canopy_factor")
   
    # clip the resulting hsi raster to the area of interest
    desc = arcpy.Describe(aoi).extent
    extent = "{} {} {} {}".format(desc.XMin,desc.YMin,desc.XMax,desc.YMax)
    arcpy.Clip_management(result[0],extent,wks+r"\coyote_hsi",aoi,"0",
                         "ClippingGeometry","MAINTAIN_EXTENT")
   
   
    
    
      










