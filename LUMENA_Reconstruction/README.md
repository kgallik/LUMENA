# Tissue reconstruction with LUMENA data

## Required material/resources:

- Python env with Python3.10 or greater
  - Use [requirements.txt](requirements.txt) to set up required dependencies
- Recommended minimum RAM: 32GB
- Downsampled images exported from QuPath with the annotation (tissue) mask and two cardinal points (East and North) as grey scale label images. Use the [batch export script](Batch_label_img_export_QUiet.groovy) modified from the [UW-LOCI's QUiet extension](https://github.com/uw-loci/qupath-extension-image-export-toolkit)
- Reference maps for reconstruction, should have distinct grey values for each region in reference map
- .csv key for the grey values corresponding to the img exports from QuPath and the grey values for the reference map regions (may need to create this by opening the reference maps in an image viewer like ImageJ or Napari, convert to grey scale and record the grey values in the key)
- A csv data frame containing the following metadata for the annotations (generated using the annotation level measurement export from QuPath with added columns):
  - `AnimalID`: Unique ID of the animal the sample originated from, used to organize exported data
  - `Genotype`: Genotype of the animal, used as an input argument to run pipeline on a specific genotype (required argument)
  - `Image`: name of the image that the exported mask comes from, expects the name here matches the name of the exported mask img from QuPath
  - `Tissue.ID`: Corresponds to the same class given to each tissue annotation in QuPath, used in for creating the pathway while searching for the mask image
    - Expected path for each mask: *'parent/path/`Tissue.ID`/`Image`.png'*
  - `Gland.side`: Used to organize images by the gland source to reconstruct the left and right glands separately
  - `MapBase`: The reference map to use for reconstruction, needs to be the same name as the image file for the reference map
  - `MappingID`: region of reference map the tissue comes from, expects this value and the reference key to be the same (i.e., a MappingID of 3 is used as the key in the grey value csv key to get the corresponding Map_Grey_Value)
  - `RegistrationLoc`: the bounding box corner to use for final placement of the tissue in the region of the reference map, must be one of the following values: top_left, top_right, top_middle, bottom_left, bottom_right, bottom_middle, right_middle, or left_middle
- A csv data frame containing the following metadata for the tiles generated with LUMENA_Tiles in QuPath. (Export the detection level measurements for all images to reconstruct into the same csv. Column headers will need to be renamed, see [renaming columns notebook](renaming_columns.ipynb) as an example for quickly renaming column headers with pandas.)
  - `Tiles_Image`: name of the image that the exported mask comes from, expects the name here matches the name of the exported mask img from QuPath
  - `Tiles_Tissue_ID`: Corresponds to the same class given to each tissue annotation in QuPath, used in for creating the pathway while searching for the mask image
  - `Tiles_Centroid_X_um`: used for calculating the transformed location after registering the tissue to the reference map
  - `Tiles_Centroid_Y_um`: used for calculating the transformed location after registering the tissue to the reference map
  - `Tiles_AnimalID`: Unique ID of the animal the sample originated from, used to organize exported data
  - `Tiles_Gland_side`: Used to organize images by the gland source to reconstruct the left and right glands separately

## Create Python env:

conda python env manager

```bash
conda create -n LUMENA python=3.10
conda activate LUMENA
#should work with new versions of python as well
```

or

venv

```bash
python3 -m venv LUMENA
source LUMENA/bin/activate
```

Install dependencies

```bash
pip install -r lumena_requirements.txt
```

## To run LUMENA reconstruction:

