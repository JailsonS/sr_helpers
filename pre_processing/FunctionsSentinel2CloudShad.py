import ee

CLD_PRB_THRESH = 60
NIR_DRK_THRESH = 0.15
CLD_PRJ_DIST = 1
BUFFER = 50

def removeShadowAndClouds2(srCollection, propCollection):
    col = getCombinedCollection(srCollection, propCollection)

    col = col.map(lambda img: addCldShdwMask(img))

    return col

def getCombinedCollection(srCollection, propCollection):


    filter = ee.Filter.equals(leftField='system:index', rightField='system:index')
    typeFilter = ee.Join.saveFirst(matchKey='s2cloudless')
    
    joined = typeFilter.apply(
        primary=srCollection,
        secondary=propCollection, 
        condition=filter
    )


    return ee.ImageCollection(joined)


def addCloudBands(img):
    # Get s2cloudless image, subset the probability band.
    cldPrb = ee.Image(img.get('s2cloudless')).select('probability')

    # Condition s2cloudless by the probability threshold valu;e.
    isCloud = cldPrb.gt(CLD_PRB_THRESH).rename('clouds')

    # Add the cloud probability layer and cloud mask as image bands.
    return img.addBands(ee.Image([cldPrb, isCloud]))


def addShadowBands(img):
    # Identify water pixels from the SCL band.
    not_water = img.select('SCL').neq(6)

    # Identify dark NIR pixels that are not water (potential cloud shadow pixels).
    SR_BAND_SCALE = 1e4
    dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels')

    # Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')))

    # Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
        .reproject(crs=img.select(0).projection(), scale=100)
        .select('distance')
        .mask()
        .rename('cloud_transform'))

    # Identify the intersection of dark pixels with cloud shadow projection.
    shadows = cld_proj.multiply(dark_pixels).rename('shadows')

    # Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))

def addCldShdwMask(img):
    # Add cloud component bands.
    img_cloud = addCloudBands(img)

    # Add cloud shadow component bands.
    img_cloud_shadow = addShadowBands(img_cloud)

    # Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
    is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0)

    # Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
    # 20 m scale is for speed, and assumes clouds don't require 10 m precision.
    is_cld_shdw = (is_cld_shdw.focalMin(2).focalMax(BUFFER*2/20)
        .reproject(
            crs=img.select([0]).projection(), 
            scale=20
        )
        .rename('cloudmask'))

    # Add the final cloud-shadow mask to the image.
    return img_cloud_shadow.addBands(is_cld_shdw)
