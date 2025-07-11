import sys
import logging
import os

from tqdm import tqdm
from shapely.wkb import loads as load_wkb
from rtree import index
from ecoshard import taskgraph
from osgeo import ogr
from osgeo import gdal

SUBWATERSHED_VECTOR_PATH = "./data/subwatershed_ca/watershed_ca.shp"
POPULATION_VECTOR_PATH = "./data/Population_ca/pop_ca_21.shp"
DEM_RASTER_PATH = "./data/dem_ca_bu.tif"
# PCUID Population


WORKING_DIR = "./workspace"
os.makedirs(WORKING_DIR, exist_ok=True)

logging.basicConfig(
    level=logging.DEBUG,
    stream=sys.stdout,
    format=(
        "%(asctime)s (%(relativeCreated)d) %(processName)s %(levelname)s "
        "%(name)s [%(funcName)s:%(lineno)d] %(message)s"
    ),
)
LOGGER = logging.getLogger(__name__)
logging.getLogger("ecoshard.taskgraph").setLevel(logging.INFO)


def build_spatial_index(layer):
    idx = index.Index()
    id_to_geom = {}
    for fid, feature in enumerate(layer):
        geometry_wkb = bytes(feature.GetGeometryRef().ExportToWkb())
        geom = load_wkb(geometry_wkb)
        idx.insert(fid, geom.bounds)
        id_to_geom[fid] = (feature.GetField("HYBAS_ID"), geom)
    return idx, id_to_geom


def build_downstream_lookup(_salt):
    ds = ogr.Open(SUBWATERSHED_VECTOR_PATH)
    layer = ds.GetLayer()

    id_to_nextdown = {}

    for feature in layer:
        hybas_id = feature.GetField("HYBAS_ID")
        next_down = feature.GetField("NEXT_DOWN")
        id_to_nextdown[hybas_id] = next_down

    # Memoization cache
    downstream_cache = {}

    # Recursive function with memoization
    def get_downstream_chain(hybas_id):
        if hybas_id in downstream_cache:
            return downstream_cache[hybas_id]

        next_down = id_to_nextdown.get(hybas_id)
        if next_down is None or next_down == 0:
            downstream_cache[hybas_id] = []
            return []

        chain = [next_down] + get_downstream_chain(next_down)
        downstream_cache[hybas_id] = chain
        return chain

    # Build full mapping
    hybas_to_downstream = {
        hybas_id: get_downstream_chain(hybas_id) for hybas_id in id_to_nextdown
    }
    return hybas_to_downstream


def generate_downstream_shapefile(hybas_to_downstream, target_hybas_id, output_path):
    ds = ogr.Open(SUBWATERSHED_VECTOR_PATH)
    layer = ds.GetLayer()
    srs = layer.GetSpatialRef()

    id_to_nextdown = {}
    id_to_feature = {}

    for feature in layer:
        hybas_id = feature.GetField("HYBAS_ID")
        next_down = feature.GetField("NEXT_DOWN")
        id_to_nextdown[hybas_id] = next_down
        id_to_feature[hybas_id] = feature.Clone()

    driver = ogr.GetDriverByName("GPKG")
    out_ds = driver.CreateDataSource(output_path)
    out_layer = out_ds.CreateLayer("downstream_chain", srs, ogr.wkbPolygon)

    # Add fields
    layer_defn = layer.GetLayerDefn()
    for i in range(layer_defn.GetFieldCount()):
        field_defn = layer_defn.GetFieldDefn(i)
        out_layer.CreateField(field_defn)

    out_defn = out_layer.GetLayerDefn()

    # Collect target and downstream features
    ids_to_write = [target_hybas_id] + hybas_to_downstream.get(target_hybas_id, [])

    for hybas_id in ids_to_write:
        original_feature = id_to_feature.get(hybas_id)
        if original_feature:
            new_feature = ogr.Feature(out_defn)
            # Copy fields
            for i in range(out_defn.GetFieldCount()):
                new_feature.SetField(
                    out_defn.GetFieldDefn(i).GetNameRef(),
                    original_feature.GetField(i),
                )
            # Copy geometry
            geom = original_feature.GetGeometryRef()
            new_feature.SetGeometry(geom.Clone())
            out_layer.CreateFeature(new_feature)
            new_feature = None  # free resources

    out_ds = None


def analyze_population_efficient(hybas_to_downstream):
    print("open population layer")
    pop_ds = ogr.Open(POPULATION_VECTOR_PATH)
    pop_layer = pop_ds.GetLayer()

    print("open subwatershed layer")
    sub_ds = ogr.Open(SUBWATERSHED_VECTOR_PATH)
    sub_layer = sub_ds.GetLayer()
    print("build spatial index")
    spatial_index, id_to_geom = build_spatial_index(sub_layer)

    hybas_to_pcuids = {}
    pcuid_to_population = {}

    total_features = pop_layer.GetFeatureCount()

    for pop_feature in tqdm(
        pop_layer, total=total_features, desc="Processing population features"
    ):
        pcuid = pop_feature.GetField("PCUID")
        population = pop_feature.GetField("Population")
        pcuid_to_population[pcuid] = int(population)

        pop_geom = load_wkb(bytes(pop_feature.GetGeometryRef().ExportToWkb()))
        candidates = spatial_index.intersection(pop_geom.bounds)

        for fid in candidates:
            hybas_id, sub_geom = id_to_geom[fid]
            if sub_geom.intersects(pop_geom):
                hybas_to_pcuids.setdefault(hybas_id, set()).add(pcuid)

    hybas_to_all_downstream_pcuids = {
        hybas_id: set.union(
            *(
                hybas_to_pcuids.get(down_id, set())
                for down_id in [hybas_id] + downstream_ids
            )
        )
        for hybas_id, downstream_ids in hybas_to_downstream.items()
    }

    hybas_to_downstream_population = {
        hybas_id: sum(pcuid_to_population[pcuid] for pcuid in pcuids)
        for hybas_id, pcuids in hybas_to_all_downstream_pcuids.items()
    }

    return hybas_to_downstream_population


def generate_population_vector(hybas_to_downstream_population, output_path):
    ds = ogr.Open(SUBWATERSHED_VECTOR_PATH)
    layer = ds.GetLayer()
    srs = layer.GetSpatialRef()

    id_to_nextdown = {}
    id_to_feature = {}

    for feature in layer:
        hybas_id = feature.GetField("HYBAS_ID")
        next_down = feature.GetField("NEXT_DOWN")
        id_to_nextdown[hybas_id] = next_down
        id_to_feature[hybas_id] = feature.Clone()
    layer.ResetReading()

    driver = ogr.GetDriverByName("GPKG")
    out_ds = driver.CreateDataSource(output_path)
    out_layer = out_ds.CreateLayer("subwatershed_population", srs, ogr.wkbPolygon)

    layer_defn = layer.GetLayerDefn()
    for i in range(layer_defn.GetFieldCount()):
        out_layer.CreateField(layer_defn.GetFieldDefn(i))
    out_layer.CreateField(ogr.FieldDefn("ds_population", ogr.OFTInteger))

    out_defn = out_layer.GetLayerDefn()

    for hybas_id, feature in id_to_feature.items():
        new_feature = ogr.Feature(out_defn)
        for i in range(layer_defn.GetFieldCount()):
            new_feature.SetField(
                out_defn.GetFieldDefn(i).GetNameRef(), feature.GetField(i)
            )
        new_feature.SetField(
            "ds_population", hybas_to_downstream_population.get(hybas_id, 0)
        )
        new_feature.SetGeometry(feature.GetGeometryRef().Clone())
        out_layer.CreateFeature(new_feature)

    out_ds = None


def rasterize_population_vector(population_vector_path, output_raster_path):
    dem_ds = gdal.Open(DEM_RASTER_PATH)
    geotransform = dem_ds.GetGeoTransform()
    projection = dem_ds.GetProjection()
    x_size = dem_ds.RasterXSize
    y_size = dem_ds.RasterYSize

    source_ds = ogr.Open(population_vector_path)
    source_layer = source_ds.GetLayer()

    options = ["ATTRIBUTE=ds_population", "COMPRESS=LZW"]
    target_ds = gdal.GetDriverByName("GTiff").Create(
        output_raster_path, x_size, y_size, 1, gdal.GDT_Float64, options
    )
    target_ds.SetGeoTransform(geotransform)
    target_ds.SetProjection(projection)

    gdal.RasterizeLayer(
        target_ds, [1], source_layer, options=["ATTRIBUTE=ds_population"]
    )

    target_ds = None


def main():
    task_graph = taskgraph.TaskGraph(WORKING_DIR, -1)
    print("build downstream lookup")
    build_downstream_task = task_graph.add_task(
        func=build_downstream_lookup,
        args=(0,),
        store_result=True,
        task_name="build downstream lookup",
    )
    hybas_to_downstream = build_downstream_task.get()

    print("analyze pop")
    analyze_pop_task = task_graph.add_task(
        func=analyze_population_efficient,
        args=(hybas_to_downstream,),
        store_result=True,
        task_name="analyze population",
    )

    print("generate downstream population vector")
    hybas_to_downstream_population = analyze_pop_task.get()
    generate_population_vector(hybas_to_downstream_population, "ds_pop.gpkg")

    rasterize_population_vector("ds_pop.gpkg", "ds_pop.tif")

    # generate_downstream_shapefile(build_downstream_task.get(), 7120311190, "ds.gpkg")
    task_graph.join()
    task_graph.close()


if __name__ == "__main__":
    main()
