#!/usr/bin/env python3 

import svg
import sys

kelly_colors = ['F3C300', '875692', 'F38400', 'A1CAF1', 'BE0032', 'C2B280', '848482', '008856', 'E68FAC', '0067A5', 'F99379', '604E97', 'F6A600', 'B3446C', 'DCD300', '882D17', '8DB600', '654522', 'E25822', '2B3D26']

def make_drawing(height=20, *args):

    # Keep track of where we are drawing
    y_pos = 0
    
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
                    x=start, y=y_pos,
                    width = end-start, height=height,
                    fill=f"#{color}"
                )
            )
            x_max = max(x_max, end)
        # Shift to a new 'row' in the figure
        y_pos += height *2
        y_max = y_pos + height

    return svg.SVG(width=x_max, height=y_max, elements=elements)

    

if __name__ == '__main__':

    l1 = [(0, 100), (340, 500)]
    l2 = [(150, 200), (340, 430)]
    t = [l1, l2]
    print(t)
    print(make_drawing(50, *t))

