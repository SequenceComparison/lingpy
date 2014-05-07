import colorsys

def colorRange(
        number,
        brightness = 300,
        ):
    """
    Function returns different colors for the given range.

    Notes
    -----

    Idea taken from
    http://stackoverflow.com/questions/876853/generating-color-ranges-in-python .
    """

    # get the hsv codes for the given range
    hsv = [(x * 1.0 / number, 0.5, 0.5) for x in range(number)]
    
    # convert hsv to rgb
    rgb = list(map(lambda x: colorsys.hsv_to_rgb(*x), hsv))
    
    # convert rgb to html
    for i in range(number):
        rgb[i] = tuple([int(rgb[i][k] * brightness) for k in range(3)])
    html = []
    for i in range(number):
        html.append("#%02x%02x%02x" % rgb[i])

    return html

def newColorRange(
        number
        ):
    
    letters = "01234567890abcdef"

    addingslow = []
    reducingslow = []
    addingfast = []
    reducingfast = []
    
    for lA in letters[::-1]:
        for lB in letters[::-1]:
            reducingslow.append(lA+lB)
    for lA in letters:
        for lB in letters:
            addingslow.append(lA+lB)
    for i in range(0,16,2):
        for j in range(0,16,2):
            addingslow.append(letters[i]+letters[j])
    for i in range(15,0,-2):
        for j in range(0,16,2):
            reducingslow.append(letters[i]+letters[j])


    colors = []
    # start with blue adding green
    for l in addingslow:
        colors.append('#00{0}ff'.format(l))
    # reduce blue
    for l in reducingfast:
        colors.append('#00ff{0}'.format(l))
    # add red
    for l in addingslow:
        colors.append('#{0}ff00'.format(l))
    # reduce green
    for l in reducingslow:
        colors.append('#ff{0}00'.format(l))

    # determine the steps for the colors
    steps = len(colors) / number

    return [colors[i] for i in range(0,len(colors),steps)]
