def read_glyph():
    from xml.dom import minidom
    glyph_path = {}
    for i in range(65,91):
        svg_file_name = 'svg/'+str(i)+'.svg'
        svg = minidom.parse(svg_file_name)
        glyph_path[chr(i)]= [path.getAttribute('d') for path in svg.getElementsByTagName('path')][-1]
    return glyph_path




