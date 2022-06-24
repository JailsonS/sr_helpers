import ee

ENDMEMBERS = [
    [119.0, 475.0, 169.0, 6250.0, 2399.0, 675.0],  # gv
    [1514.0, 1597.0, 1421.0, 3053.0, 7707.0, 1975.0],  # npv
    [1799.0, 2479.0, 3158.0, 5437.0, 7707.0, 6646.0],  # soil
    [4031.0, 8714.0, 7900.0, 8989.0, 7002.0, 6607.0]  # loud
]

def addNdvi(image: ee.image.Image) -> ee.image.Image:
    ndvi = image.normalizedDifference(['nir','red']).float()
    return image.addBands(ndvi.rename('ndvi'))


def addNdwi(image: ee.image.Image) -> ee.image.Image:
    ndwi = image.normalizedDifference(['nir','swir1']).float()
    return image.addBands(ndwi.rename('ndwi'))


def addNdfi(image: ee.image.Image) -> ee.image.Image:

    summed = image.expression('b("gv") + b("npv") + b("soil")')

    gvs = image.select("gv").divide(summed).multiply(100).byte().rename("gvs")

    npvSoil = image.expression('b("npv") + b("soil")').rename('npvSoil')

    ndfi = ee.Image.cat(gvs, npvSoil).normalizedDifference().rename('ndfi')

    # rescale NDFI from 0 to 200 \
    ndfi = ndfi.expression('byte(b("ndfi") * 100 + 100)').rename('ndfi')
    image = image.addBands(gvs)
    image = image.addBands(ndfi)

    return image


def addCsfi(image: ee.image.Image) -> ee.image.Image:
    csfi = image.expression(
        "(float(b('gv') - b('shade'))/(b('gv') + b('shade')))")

    csfi = csfi.multiply(100).add(100).byte().rename(['csfi'])

    return image.addBands(csfi)


def getFranctions(image: ee.image.Image) -> ee.image.Image:
    
    outBandNames = ['gv', 'npv', 'soil', 'cloud']
    
    fractions = ee.Image(image)\
        .select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])\
        .unmix(ENDMEMBERS)\
        .max(0)\
        .multiply(100)\
        .byte()
    
    fractions = fractions.rename(outBandNames)
    
    summed = fractions.expression('b("gv") + b("npv") + b("soil")')
    
    shade = summed.subtract(100).abs().byte().rename("shade")

    fractions = fractions.addBands(shade)
    fractions = fractions.copyProperties(image)\
        .copyProperties(image, [
            'system:time_start',
            'system:time_end',
            'system:footprint'
    ])
    
    return image.addBands(fractions)


def addEvi(image: ee.image.Image) -> ee.image.Image:
    exp = "2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)"

    evi = image.expression(exp, {
        'nir': image.select('nir'),
        'red': image.select('red'),
        'blue': image.select('blue')
    })

    evi = evi.rename("evi").float()

    return image.addBands(evi)