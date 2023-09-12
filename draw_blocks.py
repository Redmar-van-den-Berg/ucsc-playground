#!/usr/bin/env python3 

import svg
import sys

kelly_colors = ['F3C300', '875692', 'F38400', 'A1CAF1', 'BE0032', 'C2B280', '848482', '008856', 'E68FAC', '0067A5', 'F99379', '604E97', 'F6A600', 'B3446C', 'DCD300', '882D17', '8DB600', '654522', 'E25822', '2B3D26']

def make_drawing(height=20, *args):

    # Keep track of where we are drawing
    y_pos = 0
    x_pos = 100
    
    # Keep track of the maximum value we drew
    x_max = 0
    y_max = height

    # Here we keep track of the elements of the drawing
    elements = list()

    for i, blocks in enumerate(args):
        color = kelly_colors[i%len(kelly_colors)]
        for start, end in blocks:
            elements.append(
                svg.Rect(
                    x=x_pos + start, y=y_pos,
                    width = end-start, height=height,
                    fill=f"#{color}"
                )
            )
            x_max = max(x_max, end + x_pos)
        # Shift to a new 'row' in the figure
        y_pos += height *2
        x_pos = 100
        y_max = y_pos + height

    return svg.SVG(width=x_max + 100, height=y_max, elements=elements)

def draw_regions(genomic_regions, offset):
    """Draw a dictionary of regions"""

    height = 10
    scale = 0.01
    elements = list()

    x_pos = 100
    y_pos = 0

    x_max = x_pos
    y_max = y_pos

    color_counter = 0
    for track, regions in genomic_regions.items():
        color = kelly_colors[color_counter%len(kelly_colors)]
        # First, we write the name of the track

        # Next, we draw Regions
        for region in regions:
            # First the blocks
            elements.append(
                svg.Rect(
                    x=(x_pos + region.start - offset) * scale, y=y_pos,
                    width = region.size * scale, height=height,
                    fill=f"#{color}"
                )
            )

            # Then below, the name of every region
            elements.append(
                svg.Text(
                    x=(x_pos + region.start - offset + region.size*0.5)*scale,
                    y=y_pos + 0.85* height,
                    text=region.name
                )
            )
            x_max = max(x_max, region.start - offset + region.size)

        y_pos += height *2
        x_pos = 100
        y_max = y_pos + height

        color_counter += 1

    return svg.SVG(width=x_max*scale + 100, height=y_max, elements=elements)

    

if __name__ == '__main__':

    l1 = [(0, 100), (340, 500)]
    l2 = [(150, 200), (340, 430)]
    t = [l1, l2]
    print(t)
    print(make_drawing(50, *t))

