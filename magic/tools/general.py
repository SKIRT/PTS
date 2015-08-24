

## This function ...
def split_xyz(data, mask=None):

    # Initialize lists to contain the x, y and z values
    x_values = []
    y_values = []
    z_values = []

    # Loop over all x and y values
    for x in range(data.shape[1]):
        for y in range(data.shape[0]):

            # If no mask is specified or the pixel is not masked, add the coordinates and value to the appropriate lists
            if mask is None or not mask[y,x]:

                x_values.append(x)
                y_values.append(y)
                z_values.append(data[y,x])

    # Return the lists of x, y and z values
    return x_values, y_values, z_values

## This function ...
def average_stddev(model):

    return 0.5*(model.x_stddev.value + model.y_stddev.value)