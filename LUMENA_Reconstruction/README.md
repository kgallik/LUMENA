# Tissue reconstruction with LUMENA data

Requirements:

- Python env with Python3.10 or greater
  - Use [requirements.txt](requirements.txt) to set up required dependencies
- Recommended minimum RAM: 32GB
- Downsampled images exported from QuPath with the annotation (tissue) mask and two cardinal points (East and North) as greyscale label images
- Reference maps for reconstruction, should be distinct grey values for each region in reference map
- .csv key for the grey values corresponding to the img exports from QuPath and the grey values for the reference map regions (may need to create this by opening the reference maps in an image viewer like ImageJ or Napari, convert to greyscale and record the grey values in the key)
- A csv dataframe containing the following metadata for the annotations:
  - `Image`: name of the image that the exported mask comes from, expects the name here matches the name of the exported mask img from QuPath
  - `Tissue.ID`: 