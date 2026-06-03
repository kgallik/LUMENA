import numpy as np
from PIL import Image

def detect_cardinal_point(arr, grey_value):
    """
    Detect the centroid of a cardinal point marker by its grey value.
    
    Returns (x, y) in pixel coordinates, or None if not found.
    """
    
    # Create binary mask for this grey value
    mask = (arr * (arr == grey_value)).astype(int)
        
    if not mask.any():
        print(f"  WARNING: No pixels found for grey value {grey_value} in {image_path}")
        return None
    
    # Label connected components and get centroid
    labeled_img = sk.measure.label(mask)
    props = sk.measure.regionprops(labeled_img)
    cy,cx = props[0].centroid
    return np.array([cy, cx])  # return as (y,x)

def get_cardinal_points(img, grey_values_dict):
    """
    Extract both cardinal points from an image.
    Returns dict with 'north' and 'east' as (y,x) arrays.
    """
    points = {}
    for direction in ["north", "east"]:
        grey = grey_values_dict[direction]
        pt = detect_cardinal_point(img, grey)
        if pt is None:
            raise ValueError(f"Could not find '{direction}' point in {image_path}")
        points[direction] = pt
        print(f"  {direction}: pixel (y,x) ({pt[0]:.1f}, {pt[1]:.1f})")
    return points

def compute_transform(mask_points,moving_image):
    """
    Determine the flip and/or rotation needed to orient the moving image
    to match the map, using the N and E cardinal point locations.

    Detection logic (image coordinates, y increases downward):
        North marker should be visually "above" East  → north_y < east_y
        East  marker should be visually "right" of North → east_x > north_x

    Flip cases:
        north_is_up  and east_is_right  → none
        north_is_up  and east_is_left → horizontal flip
        north_is_right and east_is_up  → -90 rotation
        north_is_right and east_is_down → 90 rotation
        north_is_left and east_is_up → Vertical flip + 90 rotation
        north_is_left and east_is_down → Vertical flip + -90 rotation
        north_is_down and east_is_left → vertical + horizontal flip
        north_is_down and east_is_right → vertical flip

    Parameters
    ----------
    mask_points  : dict with 'north' and 'east' as (x, y) pixel arrays — moving image
    moving_image: image array of moving image
    Returns
    -------
    flip matrix
    """
    # --- Detect required flip from cardinal point geometry ---
    # Compare pixel coords directly for orientation — spacing does not
    # affect which direction is "up" or "right"
    #get moving image relative cardinal coordinates:
    north = (0, moving_image.shape[1]//2)
    south = (moving_image.shape[0]-1, moving_image.shape[1]//2)
    east = (moving_image.shape[0]//2, moving_image.shape[1]-1)
    west = (moving_image.shape[0]//2, 0)
    cardinals = {"north": north, "south": south, "east": east, "west": west}

    #which cardinal direction are the two reference points closest to?
    closest_north = min(cardinals, 
                  key=lambda d: np.hypot(mask_points['north'][0] - cardinals[d][0], 
                    mask_points['north'][1] - cardinals[d][1]))
    closest_east = min(cardinals, 
                  key=lambda d: np.hypot(mask_points['east'][0] - cardinals[d][0], 
                    mask_points['east'][1] - cardinals[d][1]))

    if closest_north == 'north' and closest_east == 'east':
        flip = None
        rotation = None
    elif closest_north == 'north' and closest_east == 'west':
        flip = 'horizontal'
        rotation = None
    elif closest_north == 'west' and closest_east == 'north':
        flip = None
        rotation = 'clockwise'
    elif closest_north == 'east' and closest_east == 'north':
        flip = 'vertical'
        rotation = 'counter-clockwise'
    elif closest_north == 'west' and closest_east == 'south':
        flip = 'vertical'
        rotation = 'clockwise'
    elif closest_north == 'east' and closest_east == 'south':
        flip = None
        rotation = 'counter-clockwise'
    elif closest_north == 'south' and closest_east == 'west':
        flip = 'both'
        rotation = None
    elif closest_north == 'south' and closest_east == 'east':
        flip = 'vertical'
        rotation = None
    else:
        print('Orientation does not match conditions, check image')

    return flip, rotation