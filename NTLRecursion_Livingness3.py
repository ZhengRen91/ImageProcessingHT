import os
from matplotlib import pyplot as plt
os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd
from statistics import mean
from shapely.geometry import shape
import fiona
import datetime
import multiprocessing
from shapely import speedups
import pandas as pd
import rasterio
from rasterio import features
import htb, htb2, htb3
from pyproj import CRS

import powerlaw
import matplotlib.pyplot as plt

starttime = datetime.datetime.now()

def load_shp(filename):
    geomlist = []
    codelist = []
    idlist = []
    with fiona.open(filename, 'r') as shp:
        for feature in shp:
            geom = shape(feature['geometry'])
            prop = feature['properties']
            geomlist.append(geom)
            codelist.append(prop['gridcode'])
            # codelist.append(prop['GRIDCODE'])
            # idlist.append(prop['ID'])
            idlist.append(prop['Id'])
    print(len(geomlist))
    shpdf = gpd.GeoDataFrame(geometry=geomlist)
    shpdf['id'] = idlist
    shpdf['gridcode'] = codelist
    # print(shpdf)
    return shpdf


def load_shp2(filename):
    geomlist = []
    codelist = []
    idlist = []
    with fiona.open(filename, 'r') as shp:
        for feature in shp:
            geom = shape(feature['geometry'])
            prop = feature['properties']
            geomlist.append(geom)
            codelist.append(prop['GRIDCODE'])
            idlist.append(prop['ID'])
    print(len(geomlist))
    shpdf = gpd.GeoDataFrame(geometry=geomlist)
    shpdf['id'] = idlist
    shpdf['gridcode'] = codelist
    print('In total: '+str(len(shpdf))+' grids loaded')
    return shpdf


def cluster_polygon(polygondf, itertimes):
    # checkvalid
    for p in polygondf.geometry:
        if p.is_empty is True:
            print (p)
    clusters = gpd.GeoDataFrame(geometry=list(polygondf.unary_union.geoms))
    clusters['clusterid'+str(itertimes)] = clusters.index+1
    clusterresult = gpd.tools.sjoin(polygondf, clusters, how='left').drop(columns='index_right')
    # clusterresult = gpd.tools.sjoin(polygondf, clusters, how='left').drop('index_right', 1)

    # print(clusterresult)
    grouped = clusterresult.groupby(['clusterid' + str(itertimes)])
    # clusterresult.plot(column = 'clusterid'+str(itertimes), cmap = "rainbow")
    # plt.show()
    return grouped


def iterate_clustering(groupeddf):
    # load first iteration
    featurelist = []
    for grids in groupeddf:
        featurelist.append(grids[1])
    # print(featurelist[0])
    # print(len(featurelist))
    lit_pixels = []
    # start division
    for subs in featurelist:
        cutoff = mean(subs['gridcode'])
        # print(cutoff)
        # iterate df itertuples much faster than iterrows
        for row in subs.itertuples(index=False):
            if row.gridcode > cutoff:
                lit_pixels.append(row)
    lit_pixelsdf = gpd.GeoDataFrame(lit_pixels)
    # print(lit_pixelsdf)
    return lit_pixelsdf


def iterate_division(pixeldf):
    lit_pixels = []
    # start division
    cutoff = mean(pixeldf['gridcode'])
    # print(cutoff)
    # iterate df itertuples much faster than iterrows
    for row in pixeldf.itertuples(index=False):
        if row.gridcode > cutoff:
            lit_pixels.append(row)
    lit_pixelsdf = gpd.GeoDataFrame(lit_pixels)
    # print(lit_pixelsdf)
    return lit_pixelsdf


def group_decomposable(groupeddf):
    featurelist = []
    pixellist = []
    for grids in groupeddf:
        featurelist.append(grids[1])
    # print(len(featurelist))
    for subs in featurelist:
        for row in subs.itertuples(index=False):
            pixellist.append(row)
    pixeldf = gpd.GeoDataFrame(pixellist)
    # summary = pixeldf.groupby('clusterid1').agg('clusterid2')
    # print(summary.describe()['count'])
    return pixeldf


def createGraph(resultdf):
    print('create graph file from the results dataframe')


def removeidentical(resultdf):
    # There are possible duplicates after conducting the unary_union
    newresultdf = []
    for rdf in resultdf:
        # print('original number of records: '+str(len(rdf)))
        # rdf = rdf.drop_duplicates(subset=['id'])
        rdf = rdf.drop_duplicates(subset=['id'])
        # print('new number of records: '+str(len(rdf)))
        newresultdf.append(rdf)
    # print(len(newresultdf))

    return newresultdf

# start running whole datasets
def initialdivision(inputshpfile):
    gdf = gpd.read_file(inputshpfile)
    print(len(gdf))
    gdfsort = gdf.sort_values(by=['UID'])

    NCpixellist = []
    gdfgroup = gdfsort.groupby('UID')
    for nc in gdfgroup:
        sdf = nc[1]
        # print(sdf)
        NCpixellist.append(sdf)
    print(len(NCpixellist))
    return NCpixellist


def testdataset():
    shpfile = r'C:\zhgren\LivingCitiesGlobal_NTL\data\Stockholm_ntl_nc_polygon1.shp'
    shpdf = load_shp(shpfile)
    resultdf = []
    grouped = cluster_polygon(shpdf, 1)
    resultdf.append(group_decomposable(grouped))
    lit = iterate_clustering(grouped)
    grouped = cluster_polygon(lit, 2)
    resultdf.append(group_decomposable(grouped))
    print(grouped.describe()['clusterid1']['count'].max())
    for i in range(3, 100):
        if grouped.describe()['clusterid1']['count'].max() > 3:
            lit = iterate_clustering(grouped)
            grouped = cluster_polygon(lit, i)
            resultdf.append(group_decomposable(grouped))
        else:
            break
    print(resultdf)
    print(len(resultdf))
    newresultdf = removeidentical(resultdf)
    level = 1
    # for df in newresultdf:
    #     df.to_file(r'C:\zhgren\LivingCitiesGlobal_NTL\testdata\Stolevel' + str(level) + '.shp',
    #                driver="ESRI Shapefile")
    #     level += 1
    # print('export successfully!')
    fatherid = 0
    chilid = 1
    treelist = [fatherid, chilid]
    treelistlist = []
    uniquelist = []
    for df in newresultdf:
        print(df['clusterid'+str(level)].unique)
        uniquelist.append(df['clusterid'+str(level)].unique())
        level += 1


def compute_vitality(wholedataset,recsv,scsv):
    nclist = initialdivision(wholedataset)
    resultlistlist = []
    resultlist = []
    numberofsubs = []
    sublistlist = []
    vitalitylist = []
    # Iterate nc
    for ncpixel in nclist:
        resultlist = []
        numberofsubs = []
        litdf = iterate_division(ncpixel)
        # level = 0
        # sub1 = cluster_polygon(litdf, level)
        # clusterlist = []
        # for s in sub1:
        #     # print(s[1]['geometry'])
        #     for p in s[1].itertuples(index=False):
        #         clusterlist.append(p)
        # final = gpd.GeoDataFrame(clusterlist)
        if (litdf.unary_union.geom_type == 'Polygon'):
            resultlist.append(litdf)
            print(resultlist)
            resultlistlist.append(resultlist)
            sublistlist.append([1])
            vitalitylist.append(1)
            print(1)

        else:
            grouped = cluster_polygon(litdf, 1)
            numberofsubs.append(grouped.ngroups)
            # print(grouped.ngroups)
            resultlist.append(group_decomposable(grouped))
            # print(resultlist[0])
            lit = iterate_clustering(grouped)
            grouped = cluster_polygon(lit, 2)
            numberofsubs.append(grouped.ngroups)
            # print(grouped.ngroups)
            resultlist.append(group_decomposable(grouped))

            # print(grouped.describe()['clusterid1']['count'].max())
            for i in range(3, 100):
                if (grouped.describe()['clusterid1']['count'].max() > 3):
                    lit = iterate_clustering(grouped)
                    # if there no lit pixels stop
                    if (lit.empty):
                        break
                    if (len(lit) <= 3):
                        break
                    if (lit.unary_union.geom_type == 'Polygon'):
                        break
                    else:
                        # print(lit.unary_union.geom_type)
                        grouped = cluster_polygon(lit, i)
                        resultlist.append(group_decomposable(grouped))
                        numberofsubs.append(grouped.ngroups)
                    # print(grouped.ngroups)
                else:
                    break

            newresultdf = removeidentical(resultlist)
            resultlistlist.append(newresultdf)

            vitality = 0
            recursion = 1
            for s in numberofsubs:
                vitality += s * recursion
                recursion += 1
            vitalitylist.append(vitality)
            print(vitality)
            sublistlist.append(numberofsubs)

    print(len(resultlistlist))

    file = pd.DataFrame(vitalitylist)
    file.to_csv(recsv)
    sfile = pd.DataFrame(sublistlist)
    sfile.to_csv(scsv)


def generateTree():
    shpfile = r'C:\zhgren\LivingCitiesGlobal_NTL\data\Stockholm_ntl_nc_polygon1.shp'
    shpdf = load_shp(shpfile)
    resultdf = []
    grouped = cluster_polygon(shpdf, 1)
    resultdf.append(group_decomposable(grouped))
    lit = iterate_clustering(grouped)
    grouped = cluster_polygon(lit, 2)
    resultdf.append(group_decomposable(grouped))
    print(grouped.describe()['clusterid1']['count'].max())
    for i in range(3, 100):
        if grouped.describe()['clusterid1']['count'].max() > 3:
            lit = iterate_clustering(grouped)
            grouped = cluster_polygon(lit, i)
            resultdf.append(group_decomposable(grouped))
        else:
            break
    print(resultdf)
    print(len(resultdf))
    newresultdf = removeidentical(resultdf)
    level = 1
    fatherid = 0
    chilid = 1
    treelist = [fatherid, chilid]
    treelistlist = []
    uniquelist = []

    # for df in newresultdf:
    #     childids = df['clusterid' + str(level)].tolist()
    #     fartherid = fatherid + len(uniquelist)
    #     level += 1


def Initial_raster(InputRaster):
    print('start reading image:...')
    with rasterio.open(InputRaster) as src:
        print(src.meta)
        r = src.read()
        rwidth = src.width
        rheight = src.height
        print (r[0][0][0], rwidth, rheight)
        pixellist = []
        for h in range(rheight):
            for w in range(rwidth):
                # nodata value is 65535
                if r[0][h][w] != 65535.0:
                    pixellist.append(int(r[0][h][w]))
        print('number of pixels without nodata: '+str(len(pixellist)))
        htresult = htb2.htb(pixellist)
        print(htresult)
        htindex = len(htresult)+1
        # using second mean value as cutoff as initial division
        cutoff = htresult[1]
        # pixels to polygons
        featurelist = []
        valuelist = []
        shapes = rasterio.features.shapes(r[0], transform=src.transform)
        for poly, value in shapes:
            if value != 65535:
                if value >= cutoff:
                    featurelist.append(shape(poly))
                    valuelist.append(value)
        print('number of polygons: '+str(len(featurelist)))
        # polygons to geodataframe
        rsgdf = gpd.GeoDataFrame({'geometry': featurelist, 'gridcode': valuelist})
        # gdf.plot(column='value', cmap='rainbow', legend=True)
        # plt.show()
        return rsgdf

def Initial_raster2(InputRaster):
    print('start reading image:...')
    with rasterio.open(InputRaster) as src:
        print(src.meta)
        r = src.read()
        rwidth = src.width
        rheight = src.height
        print (r[0][0][0], rwidth, rheight)
        pixellist = []
        for h in range(rheight):
            for w in range(rwidth):
                # nodata value is 65535
                if r[0][h][w] !=0:
                    pixellist.append(int(r[0][h][w]))
        print('number of pixels without nodata: '+str(len(pixellist)))
        htresult = htb2.htb(pixellist)
        print(htresult)
        htindex = len(htresult)+1
        # using second mean value as cutoff as initial division
        cutoff = htresult[0]
        # pixels to polygons
        featurelist = []
        valuelist = []
        shapes = rasterio.features.shapes(r[0], transform=src.transform)
        for poly, value in shapes:
            if value != 0:
                if value >= cutoff:
                    featurelist.append(shape(poly))
                    valuelist.append(value)
        print('number of polygons: '+str(len(featurelist)))
        # polygons to geodataframe
        rsgdf = gpd.GeoDataFrame({'geometry': featurelist, 'gridcode': valuelist})
        # gdf.plot(column='value', cmap='rainbow', legend=True)
        # plt.show()
        return rsgdf

def Initial_rasterall(InputRaster):
    print('start reading image:...')
    with rasterio.open(InputRaster) as src:
        print(src.meta)
        r = src.read()
        rwidth = src.width
        rheight = src.height
        print (r[0][0][0], rwidth, rheight)
        pixellist = []
        for h in range(rheight):
            for w in range(rwidth):
                # nodata value is 65535
                if r[0][h][w] != 0:
                    pixellist.append(int(r[0][h][w]))
        print('number of pixels without nodata: '+str(len(pixellist)))
        # htresult = htb2.htb(pixellist)
        # print(htresult)
        # htindex = len(htresult)+1
        # # using second mean value as cutoff as initial division
        # cutoff = htresult[1]
        # pixels to polygons
        featurelist = []
        valuelist = []
        shapes = rasterio.features.shapes(r[0], transform=src.transform)
        for poly, value in shapes:
            if value != 0:
                featurelist.append(shape(poly))
                valuelist.append(value)
        print('number of polygons: '+str(len(featurelist)))
        # polygons to geodataframe
        rsgdf = gpd.GeoDataFrame({'geometry': featurelist, 'gridcode': valuelist})
        # gdf.plot(column='value', cmap='rainbow', legend=True)
        # plt.show()
        return rsgdf


def readlargeraster(InputRaster):
    print('start reading image:...')
    with rasterio.open(InputRaster) as src:
        print(src.meta)
        r = src.read()
        rwidth = src.width
        rheight = src.height
        print (r[0][0][0], rwidth, rheight)
        pixellist = []
        for h in range(rheight):
            for w in range(rwidth):
                # nodata value is 65535
                if r[0][h][w] != 65535.0:
                    pixellist.append(int(r[0][h][w]))
        print('number of pixels without nodata: '+str(len(pixellist)))
        # pixels to polygons
        featurelist = []
        valuelist = []
        shapes = rasterio.features.shapes(r[0], transform=src.transform)
        for poly, value in shapes:
            if value != 65535:
                featurelist.append(shape(poly))
                valuelist.append(value)
        print('number of polygons: '+str(len(featurelist)))
        # polygons to geodataframe
        rsgdf = gpd.GeoDataFrame({'geometry': featurelist, 'gridcode': valuelist})
        # gdf.plot(column='value', cmap='rainbow', legend=True)
        # plt.show()
        return rsgdf


def clusteringRPolygons(rasterpolyon):
    print('start clustering polygons:...')
    # grouping pixelpolygons into NC clusters
    grouped = cluster_polygon(rasterpolyon, 1)
    print('number of clusters: '+str(grouped.ngroups))
    return grouped


def decompasbleRPolygons(grouped):
    # decomposable polygons
    featurelist = []
    decomplsit = []
    valuelist = []
    subslist = []
    for grids in grouped:
        subslist.append(grids[1])
        if len(grids[1]) > 2:
            featurelist.append(grids[1])
    # print(len(featurelist))
    for feature in featurelist:
        values = feature['gridcode'].tolist()
        valuelist.append(values)
        htresult = htb2.htb(values)
        htindex = len(htresult) + 1
        if htindex > 2:
            decomplsit.append(feature)

    decompNum = len(decomplsit)
    SubNum = len(subslist)

    if (len(featurelist) != 0):
        flat_list = []
        for sublist in valuelist:
            for item in sublist:
                flat_list.append(item)
        H = len(htb2.htb(flat_list)) + 1
    else: H = 0

    # print('number of decomposable polygons: '+str(decompNum))
    return decomplsit, decompNum, H, SubNum


def decompasbleRPolygons2(grouped):
    # decomposable polygons
    dlist = []
    valuelist = []
    subslist = []
    pixellist = []

    for grids in grouped:
        subslist.append(grids[1])
    for subs in subslist:
        for row in subs.itertuples(index=False):
            pixellist.append(row)
        values = subs['gridcode'].tolist()
        valuelist.append(values)
        htresult = htb2.htb(values)
        htindex = len(htresult) + 1
        print(htindex)
        if htindex >= 3:
            dlist.append(subs)

    pixeldf = gpd.GeoDataFrame(pixellist)
    decompNum = len(dlist)
    SubNum = len(subslist)

    if (len(subslist) != 0):
        flat_list = []
        for vs in valuelist:
            for item in vs:
                flat_list.append(item)
        H = len(htb2.htb(flat_list)) + 1
    else: H = 0

    return dlist, decompNum, H, SubNum, pixeldf


def decompasbleRPolygons3(grouped):
    # decomposable polygons
    dlist = []
    valuelist = []
    subslist = []
    pixellist = []

    for grids in grouped:
        subslist.append(grids[1])
    for subs in subslist:
        for row in subs.itertuples(index=False):
            pixellist.append(row)
        values = subs['gridcode'].tolist()
        valuelist.append(values)
        htresult = htb2.htb(values)
        htindex = len(htresult)+1
        if htindex >= 3:
            dlist.append(subs)

    pixeldf = gpd.GeoDataFrame(pixellist)
    decompNum = len(dlist)
    SubNum = len(subslist)

    if (len(subslist) != 0):
        flat_list = []
        for vs in valuelist:
            for item in vs:
                flat_list.append(item)
        H = len(htb2.htb(flat_list))+1
    else: H = 0

    return dlist, decompNum, H, SubNum, pixeldf

def get_litall(decomplist):
    lit = []
    for decomp in decomplist:
        valuelist = []
        # using mean value as cutoff
        cutoff = mean(decomp['gridcode'])
        for row in decomp.itertuples(index=False):
            if row.gridcode > cutoff:
                lit.append(row)
                valuelist.append(row.gridcode)
        # if len(htb2.htb(valuelist)) < 2:
        #     # print('none decomposable reached')
        #     return None
    litdf = gpd.GeoDataFrame(lit)
    return litdf


def get_litone(decomp):
    lit = []
    valuelist = []
    # using mean value as cutoff
    cutoff = mean(decomp['gridcode'])
    for row in decomp.itertuples(index=False):
        if row.gridcode > cutoff:
            lit.append(row)
            valuelist.append(row.gridcode)
    # if len(htb2.htb(valuelist)) < 2:
    #     # print('none decomposable reached')
    #     return None
    litdf = gpd.GeoDataFrame(lit)
    return litdf


def iterate_decomposing(decomplsit):
    #     start division until the none decomposable reached (depth first iteration)
    print('start iterating clustering from iteration 2:...')
    LRlist = []
    Dlist = []
    Hlist = []
    Slist = []
    Itrlist = []
    for decomp in decomplsit:

        litdf = get_litone(decomp)
        if litdf is None:
            continue
        else:
            for r in range(2, 10):
                if (litdf.unary_union.geom_type == 'Polygon'):
                    break

                group_r = cluster_polygon(litdf, r)
                dlist, D, H, S = decompasbleRPolygons(group_r)
                Lr = H * S
                print('Lr: '+str(Lr)+' D: '+str(D)+' H: '+str(H)+' S: '+str(S))
                if D == 0:
                    break
                else:
                    LRlist.append(Lr)
                    Dlist.append(D)
                    Hlist.append(H)
                    Slist.append(S)
                    litdf = get_litall(dlist)
                    if litdf is None:
                        # print(r)
                        Itrlist.append(r)
                        break

    cout = len((LRlist))
    print('number of iterations: '+str(cout))
    sumLR = sum(LRlist)
    maxR = max(Itrlist)
    print('sumLR: '+str(sumLR))
    print('maxR: '+str(maxR))
    return sumLR



def rasterDecomposing():
    rasterfile = r'C:\zhgren\LivingCitiesGlobal_NTL\testdata\World_NTL_Mollweide_reclassify_16int.tif'
    rasterpolygons = Initial_raster(rasterfile)
    firstcluster = clusteringRPolygons(rasterpolygons)
    firstdecomplist, firstD, firstH, firstS = decompasbleRPolygons(firstcluster)
    print('firstD: '+str(firstD)+' firstH: '+str(firstH)+' firstS: '+str(firstS))
    LRrest = iterate_decomposing(firstdecomplist)
    LRall = firstH * firstS + LRrest
    print('LRall: '+str(LRall))


def largerDecomposing():
    shp = r'C:\zhgren\LivingCitiesGlobal_NTL\testdata\chris\TW_ntl_pixel_chr.shp'
    rasterpolygons = load_shp(shp)
    firstcluster = clusteringRPolygons(rasterpolygons)
    firstdecomplist, firstD, firstH, firstS = decompasbleRPolygons(firstcluster)
    print('firstD: '+str(firstD)+' firstH: '+str(firstH)+' firstS: '+str(firstS))
    LRrest = iterate_decomposing(firstdecomplist)
    LRall = firstH * firstS + LRrest
    print('LRall: '+str(LRall))


def largerDecomposingShp():
    shp = r'C:\zhgren\LivingCitiesGlobal_NTL\testdata\chris\TW_ntl_pixels.shp'
    rasterpolygons = load_shp(shp)
    firstgrouped = clusteringRPolygons(rasterpolygons)
    firstdecomplist, firstD, firstH, firstS, firstPixels = decompasbleRPolygons2(firstgrouped)
    print('firstD: '+str(firstD)+' firstH: '+str(firstH)+' firstS: '+str(firstS))
    LRlist = []
    LRlist.append(firstH * firstS)

    groupedpixels = []
    groupedpixels.append(firstPixels)
    resultlistlist = []
    depthlist = []
    # resultlistlist.append(groupedpixels)

    for decomp in firstdecomplist:
        resultlist = []
        litdf = get_litone(decomp)
        print('litnumber1: '+str(len(litdf)))
        if litdf is None:
            continue
        # if len(litdf) <4 and len(litdf) >0:
        #     continue
        else:
            for r in range(2, 10):
                print('this is iteration: '+str(r))
                if (litdf.unary_union.geom_type == 'Polygon'):
                    print('none decomposable reached')
                    break
                else:
                    group_r = cluster_polygon(litdf, r)
                    dlist, D, H, S, pixels = decompasbleRPolygons2(group_r)
                    Lr = H * S
                    print('Lr: ' + str(Lr) + ' D: ' + str(D) + ' H: ' + str(H) + ' S: ' + str(S))
                    resultlist.append(pixels)
                    LRlist.append(Lr)

                    if D == 0:
                        break
                    else:
                        litdf = get_litall(dlist)
                        print('litnumber2: ' + str(len(litdf)))
            if(len(resultlist)>0):
                resultlistlist.append(resultlist)
                depth = len(resultlist)
                depthlist.append(depth)
    print('number of iterations: '+str(len(resultlistlist)))
    print('max depth: '+str(max(depthlist)))

    firstPixels.to_file(r'C:\zhgren\LivingCitiesGlobal_NTL\data\world_mask_result\NCs\World0.shp',
                        driver="ESRI Shapefile")
    print('first level export successfully!')

    nclist = []
    for i in range(0, max(depthlist)):
        print(i)
        singlevelnc = []
        for nc in resultlistlist:
            try:
                singlevelnc.append(nc[i])
            except:
                continue
        nclist.append(singlevelnc)

    print('nc list length:'+str(len(nclist)))
    # print(nclist[0])

    for level in range(0, max(depthlist)):
        df = pd.concat(nclist[level])
        rdf = df.drop_duplicates(subset=['id'])
        rdf.to_file(r'C:\zhgren\LivingCitiesGlobal_NTL\data\world_mask_result\NCs\World' + str(level+1) + '.shp',
                    driver="ESRI Shapefile")

        print('export successfully!')

    # file = pd.DataFrame(LRlist)
    # file.to_csv(r'C:\zhgren\LivingCitiesGlobal_NTL\testdata\chris\TW_NCs\TW_LR.csv', index=False, header=False)


def largerDecomposingShp_global():
    shp = r'C:\zhgren\LivingCitiesGlobal_NTL\data\world_mask_result\world_part2a.shp'
    rasterpolygons = load_shp(shp)

    print('rasterpolygons: '+str(len(rasterpolygons)))

    firstgrouped = clusteringRPolygons(rasterpolygons)
    firstdecomplist, firstD, firstH, firstS, firstPixels = decompasbleRPolygons2(firstgrouped)
    print('firstD: '+str(firstD)+' firstH: '+str(firstH)+' firstS: '+str(firstS))
    LRlist = []
    LRlist.append(firstH * firstS)

    groupedpixels = []
    groupedpixels.append(firstPixels)
    resultlistlist = []
    depthlist = []
    # resultlistlist.append(groupedpixels)

    for decomp in firstdecomplist:
        resultlist = []
        litdf = get_litone(decomp)
        print('litnumber1: '+str(len(litdf)))
        if litdf is None:
            continue
        # if len(litdf) <4 and len(litdf) >0:
        #     continue
        else:
            for r in range(2, 10):
                print('this is iteration: '+str(r))
                if (litdf.unary_union.geom_type == 'Polygon'):
                    print('none decomposable reached')
                    break
                else:
                    group_r = cluster_polygon(litdf, r)
                    dlist, D, H, S, pixels = decompasbleRPolygons2(group_r)
                    Lr = H * S
                    print('Lr: ' + str(Lr) + ' D: ' + str(D) + ' H: ' + str(H) + ' S: ' + str(S))
                    resultlist.append(pixels)
                    LRlist.append(Lr)

                    if D == 0:
                        break
                    else:
                        litdf = get_litall(dlist)
                        print('litnumber2: ' + str(len(litdf)))

            if(len(resultlist)>0):
                resultlistlist.append(resultlist)
                depth = len(resultlist)
                depthlist.append(depth)

    print('number of iterations: '+str(len(resultlistlist)))
    print('max depth: '+str(max(depthlist)))

    firstPixels.to_file(r'C:\zhgren\LivingCitiesGlobal_NTL\data\world_mask_result\NCs2\World_part2a_0.shp',
                        driver="ESRI Shapefile")
    print('first level export successfully!')

    nclist = []
    for i in range(0, max(depthlist)):
        print(i)
        singlevelnc = []
        for nc in resultlistlist:
            try:
                singlevelnc.append(nc[i])
            except:
                continue
        nclist.append(singlevelnc)

    print('nc list length:'+str(len(nclist)))
    # print(nclist[0])

    for level in range(0, max(depthlist)):
        df = pd.concat(nclist[level])
        rdf = df.drop_duplicates(subset=['id'])
        rdf.to_file(r'C:\zhgren\LivingCitiesGlobal_NTL\data\world_mask_result\NCs2a\World_part2a_' + str(level+1) + '.shp',
                    driver="ESRI Shapefile")

        print('export successfully!')

    # file = pd.DataFrame(LRlist)
    # file.to_csv(r'C:\zhgren\LivingCitiesGlobal_NTL\testdata\chris\TW_NCs\TW_LR.csv', index=False, header=False)


def fromshp2tree():
    shpfolder = r'C:\zhgren\LivingCitiesGlobal_NTL\testdata\chris\TW_NCs'
    shpfiles = os.listdir(shpfolder)
    shpfiles.sort()
    treelist = []
    parentlabel = []
    childlabel = []
    for file in shpfiles:
        if file.endswith(".shp"):
            # print(file)
            fullpath = os.path.join(shpfolder, file)
            print(fullpath)
            iter = int(file.split('.')[0][2]) + 1
            childcol = 'clusterid' + str(iter)
            print(childcol)
            shpdf = gpd.read_file(fullpath)
            clusternum = shpdf.groupby(childcol).ngroups
            print('clusternum: '+str(clusternum))
            for i in range(0, clusternum):
                if childcol == 'clusterid1':
                    parentlabel.append('root')
                    childid = shpdf.groupby(childcol)[childcol].mean().iloc[i]
                    clabel = str(iter)+'_'+str(int(childid))
                    print(clabel)
                    childlabel.append(clabel)
                else:
                    parentcol = 'clusterid' + str(iter-1)
                    plabel = str(iter-1)+'_'+str(int(shpdf.groupby(childcol)[parentcol].mean().iloc[i]))
                    print(plabel)
                    parentlabel.append(plabel)

                    childid = shpdf.groupby(childcol)[childcol].mean().iloc[i]
                    clabel = str(iter) + '_' + str(int(childid))
                    print(clabel)
                    childlabel.append(clabel)

    print(len(parentlabel))
    print(len(childlabel))

    for i in range(0, len(parentlabel)):
        treelist.append([parentlabel[i], childlabel[i]])

    with open(r"C:\zhgren\LivingCitiesGlobal_NTL\testdata\chris\TreeList\TW_tree_b2.txt", "w") as output:
        for row in treelist:
            output.write(str(row[0]) +',' +row[1] + '\n')


def alpha_cal():
    shpfolder = r'C:\zhgren\LivingCitiesGlobal_NTL\testdata\chris\TW_NCs'
    shpfiles = os.listdir(shpfolder)
    shpfiles.sort()
    arealist = []
    for file in shpfiles:

        if file.endswith(".shp"):

            # print(file)
            fullpath = os.path.join(shpfolder, file)
            print(fullpath)
            col = 'clusterid'+str(int(file.split('.')[0][2]) + 1)
            print(col)
            shpdf = gpd.read_file(fullpath)
            testshp = shpdf.copy()
            proj_crs = CRS.from_string('ESRI:54009')
            testshp = testshp.set_crs(proj_crs)
            print(testshp.crs)
            testshp['area'] = testshp['geometry'].area
            testshp['area'] = testshp['area'].round(2)
            gr = testshp.groupby([col])['area'].sum()
            print(gr)
            for i in range(0, len(gr)):
                arealist.append(gr.iloc[i])
            print(len(arealist))

    with open(r"C:\zhgren\LivingCitiesGlobal_NTL\testdata\chris\AreaList\TW_b2.txt", "w") as output:
        for row in arealist:
            output.write(str(row) + '\n')
    print('export successfully!')

def LivingGen_R_DN(shpfile, outputncloc):
    Hlist = []
    Slist = []
    Dlist = []
    LRlist = []
    resultlist = []
    Resultlistlist = []
    depthlist = []

    shpdf = load_shp(shpfile)
    dnlist = shpdf['gridcode'].tolist()
    print(htb2.htb(dnlist))
    h1 = len(htb2.htb(dnlist))+1
    print('H1: '+str(h1))

    firstgrouped = cluster_polygon(shpdf, 1)
    s1 = firstgrouped.ngroups
    print('S1: '+str(s1))
    Lr1 = h1 * s1
    print('LR1: ' + str(Lr1))

    d_list = []
    for group in firstgrouped:
        featuredf = group[1] # get the dataframe of each group (1,[oid, geometry, gridcode, clusterid])
        dnlist = featuredf['gridcode'].tolist()
        # arealist = featuredf['geometry'].area.tolist()
        ht = len(htb3.htb(dnlist))+1 # Use DN value as the stop condition the living image is using the H of areas.
        if(ht>=3):
            d_list.append(featuredf)
    print('Iteration 1: '+str(s1))
    print('D1: '+str(len(d_list)))

    # output the first substructures
    # pixellist = []
    # for d in d_list:
    #     for row in d.itertuples(index=False):
    #         pixellist.append(row)
    # pixeldf = gpd.GeoDataFrame(pixellist)
    # pixeldf.to_file(outputncloc, driver="ESRI Shapefile")

    # start the second iteration from the decomposable polygons
    for r in range(2, 10):
        print('this is iteration: ' + str(r))
        resultlist = []
        print('number of decomp: '+str(len(d_list)))
        if (len(d_list) == 0):
            break
        else:
            litdf = get_litall(d_list)
            if len(litdf)==0:
                print('none decomposable reached')
                break
            print('number of decomp: '+str(len(litdf)))
            if (litdf.unary_union.geom_type == 'Polygon'):
                dnlist = litdf['gridcode'].tolist()
                ht = len(htb3.htb(
                    dnlist)) + 1
                if (ht >= 3):
                    resultlist.append(litdf)
                    print('continue decompose')
                    continue
                else:
                    print('none decomposable reached')
                    break
            group_r = cluster_polygon(litdf, r)
            d_list, D, H, S, pixels = decompasbleRPolygons3(group_r)
            resultlist.append(pixels)
            Lr = H * S
            print('Lr: ' + str(Lr) + ' D: ' + str(D) + ' H: ' + str(H) + ' S: ' + str(S))
            LRlist.append(Lr)
            if D==0:
                break
            # else:
                # litdf = get_litall(d_list)

        if (len(resultlist) > 0):
            Resultlistlist.append(resultlist)


    for i in range(0, len(Resultlistlist)):
        resultlist = Resultlistlist[i]
        for j in range(0, len(resultlist)):
            result = resultlist[j]
            print(result)
            result.to_file(r'C:\zhgren\LivingCitiesGlobal_NTL\data\UKcase\londonntlgrids_S' + str(i+1) + '.shp',
                        driver="ESRI Shapefile")

    for i in range(len(LRlist)):
        LR = LRlist[i]
        print('LR: ' + str(LR))
    LR = sum(LRlist) + Lr1
    print('LRall: ' + str(LR))


def ras2shp():
    rasterfile = r'C:\zhgren\LivingCitiesGlobal_NTL\data\UKcase\londonntl.tif'
    rasterpolygons = Initial_raster2(rasterfile)
    rasterpolygons.to_file(r'C:\zhgren\LivingCitiesGlobal_NTL\data\UKcase\londonntl_grids_r0.shp',
                           driver="ESRI Shapefile")


def countingPoints(points, polygons):

    # Check the structure of the dataframes
    print(polygons.head())
    print(points.head())

    # Perform spatial join to get points within each polygon
    points_in_polygons = gpd.sjoin(points, polygons, how="left", op="within")

    # Now you have a DataFrame containing points and their corresponding polygon attributes
    print(points_in_polygons.head())

    # Group by polygon attributes and count points
    points_count = points_in_polygons.groupby('polygon_attribute_column')['point_attribute_column'].count()

    # Group by polygon attributes and calculate sum of point attributes
    points_sum = points_in_polygons.groupby('polygon_attribute_column')['point_attribute_column'].sum()

def fitpowerlaw():
    # Sample data (word frequencies)
    # reader = pd.read_table(r'C:\zhgren\LivingCitiesGlobal_NTL\data\powerlawresult\world_levels\worldlevel4.txt', sep='\t', low_memory=False,header=None)
    # freqlist=reader[0].values.tolist()
    # print(freqlist)
    # print(len(freqlist))

    shppath = r'C:\zhgren\LivingCitiesGlobal_NTL\data\world_mask_result\World_level0_H1.shp'
    shapefile = gpd.read_file(shppath)
    arealist = shapefile['geometry'].area.tolist()
    arealist.sort(reverse=True)
    print(len(arealist))

    # Fit power-law distribution using maximum likelihood estimation
    fit = powerlaw.Fit(arealist, discrete=True, estimate_discrete=True)
    R, p = fit.distribution_compare('power_law', 'lognormal')
    # Print estimated parameters
    print("Alpha (exponent of the power law):", fit.alpha)
    print("Xmin (minimum value for power-law behavior):", fit.xmin)
    print("Sigma (Standard error for power-law behavior):", fit.sigma)
    print("R: ", R)
    print("p: ", p)
    # Plot the data and power-law fit
    fig = fit.plot_ccdf(color='b', linewidth=0, marker='o', markersize=5)
    fit.power_law.plot_ccdf(color='g', linestyle='--', ax=fig)
    plt.xlabel('X')
    plt.ylabel('P(X>xmin)')
    plt.title('Power-law Distribution Fit')
    plt.tight_layout()
    plt.show()

    x_values, y_values = fit.ccdf()
    df = pd.DataFrame({'x': x_values/1000000, 'y': y_values})
    df.to_csv(r'C:\zhgren\LivingCitiesGlobal_NTL\data\powerlawresult\world_levels\worldlevel0_ccdf.csv', index=False, header=False)


def compute_L():
    level1_path = r'C:\zhgren\LivingCitiesGlobal_NTL\data\world_mask_result\World_level1_H1D.shp'
    l1 = gpd.read_file(level1_path)
    levelrest_path = r'C:\zhgren\LivingCitiesGlobal_NTL\data\world_mask_result\World_H1_level2to7_joinfather.shp'
    subs = gpd.read_file(levelrest_path)
    # pixel_path = r'C:\zhgren\LivingCitiesGlobal_NTL\data\world_mask_result\world_ntl_pt_level1.shp'
    # ntlpt = gpd.read_file(pixel_path)

    SN_list = []
    H_list = []
    L_list = []
    cityid = []

    fatherlist =[]
    childlist = []

    for i in range(0, len(l1)):
        city = l1.iloc[i]
        fatherlist.append(city)
    print(len(fatherlist))

    for i in range(0, len(subs)):
        sub = subs.iloc[i]
        childlist.append(sub)
    print(len(childlist))

    for f in fatherlist:
        cityid.append(f['uid'])
        arealist = []
        for c in childlist[:]:
            if f['uid']==(c['father_id']):
                arealist.append(c['geometry'].area)

        if arealist == []:
            # print('none decomposable reached')
            S_noneRecursion = 0
            SN_list.append(S_noneRecursion)
            H_noneRecursion = 0
            H_list.append(H_noneRecursion)
            L_noneRecursion = 0
            L_list.append(L_noneRecursion)
        else:
            S_noneRecursion = len(arealist)  # number of substructures
            SN_list.append(S_noneRecursion)

            H_noneRecursion = len(htb2.htb(arealist))
            H_list.append(H_noneRecursion)

            L_noneRecursion = H_noneRecursion * S_noneRecursion
            L_list.append(L_noneRecursion)
        print(len(cityid))


    df = pd.DataFrame({'L': L_list, 'S': SN_list, 'H': H_list, 'cityid': cityid})
    df.to_csv(r'C:\zhgren\LivingCitiesGlobal_NTL\data\world_mask_result\World_level1_H1D_L.csv', index=False, header=True)
    print('export successfully!')

def compute_SV():
    level1_path = r'C:\zhgren\LivingCitiesGlobal_NTL\data\Lcomputation\L1_top100.shp'
    l1 = gpd.read_file(level1_path)
    levelrest_path = r'C:\zhgren\LivingCitiesGlobal_NTL\data\Lcomputation\L2_7_top100.shp'
    subs = gpd.read_file(levelrest_path)

    # pixel_path = r'C:\zhgren\LivingCitiesGlobal_NTL\data\world_mask_result\world_ntl_pt_level1.shp'
    # ntlpt = gpd.read_file(pixel_path)

    cityid = []
    SV_list = []

    fatherlist =[]
    for i in range(0, len(l1)):
        city = l1.iloc[i]
        fatherlist.append(city)
    print(len(fatherlist))

    split_child_df = []
    grouped = subs.groupby('level')
    for group_id, group_df in grouped:
        split_child_df.append(group_df)
    print(len(split_child_df))

    print(split_child_df[0].head())

    for f in fatherlist:
        cityid.append(f['uid'])
        arealistlist = []

        for i in range(0, len(split_child_df)):
            arealist = []
            for j in range(0, len(split_child_df[i])):
                child = split_child_df[i].iloc[j]
                if f['uid']==(child['father_id']):
                    arealist.append(child['geometry'].area)
            if arealist == []:
                print('none decomposable reached')
            else:
                arealistlist.append(arealist)
        L = 0

        for arealist in arealistlist:
            s = len(arealist)  # number of substructures
            h = len(htb2.htb(arealist))
            l = h * s
            L = L + l

        print('SVW: '+str(L))
        SV_list.append(L)

    df = pd.DataFrame({'SV': SV_list, 'cityid': cityid})
    df.to_csv(r'C:\zhgren\LivingCitiesGlobal_NTL\data\Lcomputation\L1_top100_SV.csv', index=False, header=True)
    print('export successfully!')

def compute_SVW():
    level1_path = r'C:\zhgren\LivingCitiesGlobal_NTL\data\Lcomputation\L1_top100.shp'
    l1 = gpd.read_file(level1_path)
    levelrest_path = r'C:\zhgren\LivingCitiesGlobal_NTL\data\Lcomputation\L2_7_top100.shp'
    subs = gpd.read_file(levelrest_path)

    # pixel_path = r'C:\zhgren\LivingCitiesGlobal_NTL\data\world_mask_result\world_ntl_pt_level1.shp'
    # ntlpt = gpd.read_file(pixel_path)

    cityid = []
    SV_list = []

    fatherlist =[]
    for i in range(0, len(l1)):
        city = l1.iloc[i]
        fatherlist.append(city)
    print(len(fatherlist))

    split_child_df = []
    grouped = subs.groupby('level')
    for group_id, group_df in grouped:
        split_child_df.append(group_df)
    print(len(split_child_df))

    print(split_child_df[0].head())

    for f in fatherlist:
        cityid.append(f['uid'])
        arealistlist = []

        for i in range(0, len(split_child_df)):
            arealist = []
            for j in range(0, len(split_child_df[i])):
                child = split_child_df[i].iloc[j]
                if f['uid']==(child['father_id']):
                    arealist.append(child['geometry'].area)
            if arealist == []:
                print('none decomposable reached')
            else:
                arealistlist.append(arealist)
        L = 0
        w = 1
        for arealist in arealistlist:
            s = len(arealist)  # number of substructures
            h = len(htb2.htb(arealist))
            l = h * s * w
            L = L + l
            w = w + 1
        print('SVW: '+str(L))
        SV_list.append(L)

    df = pd.DataFrame({'SV': SV_list, 'cityid': cityid})
    df.to_csv(r'C:\zhgren\LivingCitiesGlobal_NTL\data\Lcomputation\L1_top100_SV_weighted.csv', index=False, header=True)
    print('export successfully!')


def compute_LR():
    level1_path = r'C:\zhgren\LivingCitiesGlobal_NTL\data\Lcomputation\L1_top100.shp'
    l1 = gpd.read_file(level1_path)
    levelrest_path = r'C:\zhgren\LivingCitiesGlobal_NTL\data\Lcomputation\L2_7_top100.shp'
    subs = gpd.read_file(levelrest_path)

    cityid = []
    LR_list = []

    fatherlist =[]
    for i in range(0, len(l1)):
        city = l1.iloc[i]
        fatherlist.append(city)
    print(len(fatherlist))

    sublist = []
    for j in range(0, len(subs)):
        sub = subs.iloc[j]
        sublist.append(sub)
    print(len(sublist))

    childlistlist = []
    for f in fatherlist:
        cityid.append(f['uid'])
        childlist = []
        polygonlist = []
        for c in sublist[:]:
            if f['uid']==c['father_id']:
                childlist.append(c)
                polygonlist.append(c['geometry'])
        print(len(childlist))

        child_decomplist = []

        for polygon_outer in polygonlist:
            count = 0
            inner_arealist = []
            for polygon_inner in polygonlist:
                if polygon_outer.contains(polygon_inner) and polygon_outer != polygon_inner:
                    count += 1
                    inner_arealist.append(polygon_inner.area)
            if count != 0:
                # print('number of polygon contained: '+str(count))
                child_decomplist.append(inner_arealist)

        print('number of decomposable: ' + str(len(child_decomplist)))
        LR = 0
        for decomp in child_decomplist:
            s = len(decomp)  # number of substructures
            if s == 1:
                h = 1
            else:
                h = len(htb2.htb(decomp))
            l = h * s
            LR = LR + l
        print('LR: '+str(LR))
        LR_list.append(LR)

    #     if childlist != []:
    #         childlistlist.append(childlist)
    # print(len(childlistlist))

    df = pd.DataFrame({'LR': LR_list, 'cityid': cityid})
    df.to_csv(r'C:\zhgren\LivingCitiesGlobal_NTL\data\Lcomputation\L1_top100_LR.csv', index=False, header=True)
    print('export successfully!')




if __name__ == '__main__':
    cpus = multiprocessing.cpu_count()
    print(cpus)
    speedups.enabled
    print(speedups.enabled)
    print(starttime)

    # wholedataset = r'C:\zhgren\LivingCitiesGlobal_NTL\data\World_NTL_NC_H3_pixel_Join500-885.shp'
    # re_csv = r'C:\zhgren\LivingCitiesGlobal_NTL\vital500_885.csv'
    # s_csv = r'C:\zhgren\LivingCitiesGlobal_NTL\sub500_885.csv'
    # compute_vitality(wholedataset, re_csv, s_csv)
    # rasterDecomposing()

    # largerDecomposingShp_global()
    # fromshp2tree()
    # alpha_cal()

    # ras2shp()
    # shapefielloc = r'C:\zhgren\LivingCitiesGlobal_NTL\data\UKcase\londonntl_grids_r0.shp'
    # outncloc = r'C:\zhgren\LivingCitiesGlobal_NTL\data\UKcase\londonntlgrids_d1.shp'
    # LivingGen_R_DN(shapefielloc, outncloc)

    # fitpowerlaw()

    # compute_L()
    # compute_SV()
    # compute_SVW()
    compute_LR()


    print(datetime.datetime.now() - starttime)

